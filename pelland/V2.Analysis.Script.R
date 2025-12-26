### SETUP -------------------------------------------------------------------
## load packages
pacman::p_load(readxl,tidyverse,metafor,emmeans,reshape2,
               janitor,gt,writexl,splines,performance,ggdist,ggridges,
               orchaRd,parameters,ggplot2,patchwork,bayestestR)

## read in data
df<-read_xlsx("V2.PTF.Data.xlsx")

## calculate pre-post sds from ses
df$pre.sd <- ifelse(is.na(df$pre.se), df$pre.sd, df$pre.se * sqrt(df$n))
df$post.sd <- ifelse(is.na(df$post.se), df$post.sd,df$post.se * sqrt(df$n))

## convert p to t (change scores)
df$t.value <- replmiss(df$t.value, with(df, qt(p.value/2, df=n-1,lower.tail=FALSE)))

## convert t to se (change scores)
df$delta.se <- replmiss(df$delta.se, with(df, ifelse(is.na(delta.mean), 
                                                     (post.mean - pre.mean)/t.value,delta.mean/t.value)))

## make positive
df$delta.se <- ifelse(df$delta.se < 0, df$delta.se * -1, df$delta.se)

## convert ci to se (change scores)
df$delta.se <- replmiss(df$delta.se, with(df, (delta.CI.upper - delta.CI.lower)/3.92))

## convert se to sd (change scores)
df$delta.sd <- replmiss(df$delta.sd, with(df, delta.se * sqrt(n)))

## calculate pre-post correlation coefficient for those with pre, post, and delta sds
df$ri <- (df$pre.sd^2 + df$post.sd^2 - df$delta.sd^2)/(2 * df$pre.sd * df$post.sd)

## remove values outside the range of -1 to +1 as they are likely due to
## misreporting or miscalculations in original studies
df$ri <- ifelse(between(df$ri,-1,1) == FALSE, NA, df$ri)

## then we'll convert using Fishers r to z, calculate a meta-analytic point
## estimate, and impute that across the studies with missing correlations
df.for.ri<- escalc(measure = "ZCOR", ri = ri, ni = n, data = df)

Meta_ri_str <- rma.mv(yi, vi, data=df.for.ri[df.for.ri$outcome=="Strength",],
                      slab=paste(author.year),
                      random = list(~ 1 | study/group/obs), method="REML", test="t")

Meta_ri_hyp <- rma.mv(yi, vi, data=df.for.ri[df.for.ri$outcome=="Hypertrophy",],
                      slab=paste(author.year),
                      random = list(~ 1 | study/group/obs), method="REML", test="t")

## robust variance estimates
RobuEstMeta_ri_str <- robust(Meta_ri_str, df.for.ri[df.for.ri$outcome=="Strength",]$study)
RobuEstMeta_ri_hyp <- robust(Meta_ri_hyp, df.for.ri[df.for.ri$outcome=="Hypertrophy",]$study)

z2r_str<- psych::fisherz2r(RobuEstMeta_ri_str$b[1])
z2r_hyp<- psych::fisherz2r(RobuEstMeta_ri_hyp$b[1])

df[df$outcome=="Strength",]$ri <- ifelse(is.na(df[df$outcome=="Strength",]$ri),
                                         z2r_str, df[df$outcome=="Strength",]$ri)
df[df$outcome=="Hypertrophy",]$ri <- ifelse(is.na(df[df$outcome=="Hypertrophy",]$ri),
                                            z2r_hyp, df[df$outcome=="Hypertrophy",]$ri)

## estimate change score difference SD where only pre-post data available
df$delta.sd <- replmiss(df$delta.sd, with(df, sqrt(pre.sd^2 + post.sd^2 - (2*ri*pre.sd*post.sd))))

## calculate effect sizes
## SMCR
data<-escalc(measure = "SMCR",m1i=post.mean,m2i=pre.mean,sd1i = pre.sd,ni=n,ri=ri,data = df)

## ROM
data2<-escalc(measure = "ROMC",m1i=post.mean,m2i=pre.mean,sd1i=post.sd,sd2i=pre.sd,ni=n,ri=ri,data = df)

## rename columns
data2<-data2%>%
  rename(yi.rom=yi,vi.rom=vi,ri.rom=ri)

## bind df, add weights, and load classes
data<-cbind(data,subset(data2,select=c(ri.rom,yi.rom,vi.rom)))
data$weights <- (1/sqrt(data$vi))
data$weights.rom <- (1/sqrt(data$vi.rom))

## RIR for spline models
data$spline.rir<-ifelse(data$rir.bucket=="Failure",-1,data$avg.rir)

## exponentiation observed effects for % change
data$yi.rom.exp<-(exp(data$yi.rom)-1)*100

## set knots for spline models, 0 - to distinguish between studies that
## provided a failure definition (non-volitional)
knot<-0

## subset data for different outcomes
data.str<-data%>%
  subset(outcome=="Strength")

data.hyp<-data%>%
  subset(outcome=="Hypertrophy")

## set factors
primary.data<-data%>%
  mutate(across(
    c(author.year,study,group,obs,outcome,train.status,nutrition.controlled,formal.cardio.intervention,
      other.rt.training,other.activity,progression,failure.definition,alternative.set.structure,cluster.sets,
      rest.redistribution,rest.pause,set.rep.equated,upper.lower.other,train.exercise,test.exercise,
      isometric.isokinetic.isotonic,within.between.design,train.intent,RM.max.submax),factor))

primary.data.str<-data.str%>%
  mutate(across(
    c(author.year,study,group,obs,outcome,train.status,nutrition.controlled,formal.cardio.intervention,
      other.rt.training,other.activity,progression,failure.definition,alternative.set.structure,cluster.sets,
      rest.redistribution,rest.pause,set.rep.equated,upper.lower.other,train.exercise,test.exercise,
      isometric.isokinetic.isotonic,within.between.design,train.intent,RM.max.submax),factor))

primary.data.hyp<-data.hyp%>%
  mutate(across(
    c(author.year,study,group,obs,outcome,train.status,nutrition.controlled,formal.cardio.intervention,
      other.rt.training,other.activity,progression,failure.definition,alternative.set.structure,cluster.sets,
      rest.redistribution,rest.pause,set.rep.equated,upper.lower.other,train.exercise,test.exercise,
      isometric.isokinetic.isotonic,within.between.design,train.intent,RM.max.submax),factor))

## save data frames
save(primary.data,file = "primary.data.RData")

## set means for marginal slopes
primary.str.mean.rir<-mean(primary.data.str$avg.rir,na.rm=TRUE) + c(-0.5,0.5)
primary.str.mean.spline.rir<-mean(primary.data.str$spline.rir,na.rm=TRUE) + c(-0.5,0.5)
primary.str.mean.age<-mean(primary.data.str$age,na.rm=TRUE) + c(-0.5,0.5)
primary.str.mean.sex<-mean(primary.data.str$sex.percent.male,na.rm=TRUE) + c(-0.5,0.5)
primary.str.mean.weeks<-mean(primary.data.str$weeks,na.rm=TRUE) + c(-0.5,0.5)
primary.str.mean.frequency.per.muscle<-mean(primary.data.str$frequency.per.muscle,na.rm=TRUE) + c(-0.5,0.5)
primary.str.mean.load.set<-mean(primary.data.str$load.set,na.rm=TRUE) + c(-0.5,0.5)
primary.str.mean.adj.sets.week<-mean(primary.data.str$adj.sets.week,na.rm=TRUE) + c(-0.5,0.5)
primary.str.mean.str.rir<-mean(primary.data.str$str.avg.rir,na.rm=TRUE) + c(-0.5,0.5)
primary.str.mean.frequency.per.movement<-mean(primary.data.str$frequency.per.strength.me,na.rm=TRUE) + c(-0.5,0.5)
primary.str.mean.str.specfific.load<-mean(primary.data.str$str.specific.load,na.rm=TRUE) + c(-0.5,0.5)
primary.str.mean.adj.str.sets.week<-mean(primary.data.str$adj.str.specific.sets.week,na.rm=TRUE) + c(-0.5,0.5)

primary.hyp.mean.rir<-mean(primary.data.hyp$avg.rir,na.rm=TRUE) + c(-0.5,0.5)
primary.hyp.mean.spline.rir<-mean(primary.data.hyp$spline.rir,na.rm=TRUE) + c(-0.5,0.5)
primary.hyp.mean.age<-mean(primary.data.hyp$age,na.rm=TRUE) + c(-0.5,0.5)
primary.hyp.mean.sex<-mean(primary.data.hyp$sex.percent.male,na.rm=TRUE) + c(-0.5,0.5)
primary.hyp.mean.weeks<-mean(primary.data.hyp$weeks,na.rm=TRUE) + c(-0.5,0.5)
primary.hyp.mean.frequency.per.muscle<-mean(primary.data.hyp$frequency.per.muscle,na.rm=TRUE) + c(-0.5,0.5)
primary.hyp.mean.load.set<-mean(primary.data.hyp$load.set,na.rm=TRUE) + c(-0.5,0.5)
primary.hyp.mean.adj.sets.week<-mean(primary.data.hyp$adj.sets.week,na.rm=TRUE) + c(-0.5,0.5)

## relevel set rep equated so "set" is reference level
primary.data$set.rep.equated<-fct_relevel(primary.data$set.rep.equated,c("set","rep","both"))
primary.data.str$set.rep.equated<-fct_relevel(primary.data.str$set.rep.equated,c("set","rep","both"))
primary.data.hyp$set.rep.equated<-fct_relevel(primary.data.hyp$set.rep.equated,c("set","rep","both"))

## subset data for different outcomes
data.str<-data%>%
  subset(outcome=="Strength"&avg.rir<=10)

data.hyp<-data%>%
  subset(outcome=="Hypertrophy"&avg.rir<=10)

## set factors
secondary.data<-data%>%
  mutate(across(
    c(author.year,study,group,obs,outcome,train.status,nutrition.controlled,formal.cardio.intervention,
      other.rt.training,other.activity,progression,failure.definition,alternative.set.structure,cluster.sets,
      rest.redistribution,rest.pause,set.rep.equated,upper.lower.other,train.exercise,test.exercise,
      isometric.isokinetic.isotonic,within.between.design,train.intent,RM.max.submax),factor))

secondary.data.str<-data.str%>%
  mutate(across(
    c(author.year,study,group,obs,outcome,train.status,nutrition.controlled,formal.cardio.intervention,
      other.rt.training,other.activity,progression,failure.definition,alternative.set.structure,cluster.sets,
      rest.redistribution,rest.pause,set.rep.equated,upper.lower.other,train.exercise,test.exercise,
      isometric.isokinetic.isotonic,within.between.design,train.intent,RM.max.submax),factor))

secondary.data.hyp<-data.hyp%>%
  mutate(across(
    c(author.year,study,group,obs,outcome,train.status,nutrition.controlled,formal.cardio.intervention,
      other.rt.training,other.activity,progression,failure.definition,alternative.set.structure,cluster.sets,
      rest.redistribution,rest.pause,set.rep.equated,upper.lower.other,train.exercise,test.exercise,
      isometric.isokinetic.isotonic,within.between.design,train.intent,RM.max.submax),factor))

## save data frames
save(secondary.data,file = "secondary.data.RData")

## set means for marginal slopes
secondary.str.mean.rir<-mean(secondary.data.str$avg.rir,na.rm=TRUE) + c(-0.5,0.5)
secondary.str.mean.spline.rir<-mean(secondary.data.str$spline.rir,na.rm=TRUE) + c(-0.5,0.5)
secondary.str.mean.age<-mean(secondary.data.str$age,na.rm=TRUE) + c(-0.5,0.5)
secondary.str.mean.sex<-mean(secondary.data.str$sex.percent.male,na.rm=TRUE) + c(-0.5,0.5)
secondary.str.mean.weeks<-mean(secondary.data.str$weeks,na.rm=TRUE) + c(-0.5,0.5)
secondary.str.mean.frequency.per.muscle<-mean(secondary.data.str$frequency.per.muscle,na.rm=TRUE) + c(-0.5,0.5)
secondary.str.mean.load.set<-mean(secondary.data.str$load.set,na.rm=TRUE) + c(-0.5,0.5)
secondary.str.mean.adj.sets.week<-mean(secondary.data.str$adj.sets.week,na.rm=TRUE) + c(-0.5,0.5)
secondary.str.mean.str.rir<-mean(secondary.data.str$str.avg.rir,na.rm=TRUE) + c(-0.5,0.5)
secondary.str.mean.frequency.per.movement<-mean(secondary.data.str$frequency.per.strength.me,na.rm=TRUE) + c(-0.5,0.5)
secondary.str.mean.str.specfific.load<-mean(secondary.data.str$str.specific.load,na.rm=TRUE) + c(-0.5,0.5)
secondary.str.mean.adj.str.sets.week<-mean(secondary.data.str$adj.str.specific.sets.week,na.rm=TRUE) + c(-0.5,0.5)

secondary.hyp.mean.rir<-mean(secondary.data.hyp$avg.rir,na.rm=TRUE) + c(-0.5,0.5)
secondary.hyp.mean.spline.rir<-mean(secondary.data.hyp$spline.rir,na.rm=TRUE) + c(-0.5,0.5)
secondary.hyp.mean.age<-mean(secondary.data.hyp$age,na.rm=TRUE) + c(-0.5,0.5)
secondary.hyp.mean.sex<-mean(secondary.data.hyp$sex.percent.male,na.rm=TRUE) + c(-0.5,0.5)
secondary.hyp.mean.weeks<-mean(secondary.data.hyp$weeks,na.rm=TRUE) + c(-0.5,0.5)
secondary.hyp.mean.frequency.per.muscle<-mean(secondary.data.hyp$frequency.per.muscle,na.rm=TRUE) + c(-0.5,0.5)
secondary.hyp.mean.load.set<-mean(secondary.data.hyp$load.set,na.rm=TRUE) + c(-0.5,0.5)
secondary.hyp.mean.adj.sets.week<-mean(secondary.data.hyp$adj.sets.week,na.rm=TRUE) + c(-0.5,0.5)

## relevel set rep equated so "set" is reference level
secondary.data$set.rep.equated<-relevel(secondary.data$set.rep.equated,ref = 2)
secondary.data.str$set.rep.equated<-relevel(secondary.data.str$set.rep.equated,ref = 2)
secondary.data.hyp$set.rep.equated<-relevel(secondary.data.hyp$set.rep.equated,ref = 2)

## subset data for different outcomes
data.str<-data%>%
  subset(outcome=="Strength"&str.avg.rir!="")

## rcs knots
primary.str.k=attr(rms::rcs(primary.data.str$avg.rir,4),"parms")[c(2,3)]
primary.str.b=attr(rms::rcs(primary.data.str$avg.rir,4),"parms")[c(1,4)]
primary.hyp.k=attr(rms::rcs(primary.data.hyp$avg.rir,4),"parms")[c(2,3)]
primary.hyp.b=attr(rms::rcs(primary.data.hyp$avg.rir,4),"parms")[c(1,4)]

secondary.str.k=attr(rms::rcs(secondary.data.str$avg.rir,4),"parms")[c(2,3)]
secondary.str.b=attr(rms::rcs(secondary.data.str$avg.rir,4),"parms")[c(1,4)]
secondary.hyp.k=attr(rms::rcs(secondary.data.hyp$avg.rir,4),"parms")[c(2,3)]
secondary.hyp.b=attr(rms::rcs(secondary.data.hyp$avg.rir,4),"parms")[c(1,4)]

#Overall aesthetics 
TITLE<-element_text(size = 10,face = "bold",family ="sans")
TITLES <- element_text(size = 10,face = "bold",family ="sans")
AXIS <- element_text(size = 8,family ="sans")
BACKGROUND <- element_rect(fill = "White", colour = "Black")
THEME <- theme(axis.text = AXIS, axis.title = TITLES,title = TITLE,panel.background = BACKGROUND,  
               legend.position ="none", legend.title = TITLES,aspect.ratio =, panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),strip.text = TITLES,strip.background = element_rect(fill ="gray95"))

options(scipen = 999)

expo=function(x){
  (exp(x)-1)*100
}

### PRIMARY ANALYSIS --------------------------------------------------------
# CHARACTERISTICS ---------------------------------------------------------
#Visual Descriptions of Training Interventions
primary.densitites.volume<-ggplot(data = primary.data,aes(x=adj.sets.week,color=outcome,fill=outcome))+THEME+
  theme(axis.text.y = element_text(face = "bold",angle = 90,size = 8,hjust = 0.5))+
  scale_x_continuous(breaks = c(0,4,8,12,16,20,24,28,32,36,40))+
  scale_color_manual(values = c("#0000a7","#c1272d"))+
  scale_fill_manual(values = c("#0000a7","#c1272d"))+
  geom_density_ridges(aes(y=outcome),
                      jittered_points = TRUE,
                      position = "raincloud",
                      alpha = 0.7, scale = 0.5)+
  labs(x="Direct Adjusted Sets Per Week Per Muscle Group",y="Density")

primary.densitites.load<-ggplot(data = primary.data,aes(x=load.set,color=outcome,fill=outcome))+THEME+
  theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank())+
  scale_x_continuous(breaks = c(30,40,50,60,70,80,90))+
  scale_color_manual(values = c("#0000a7","#c1272d"))+
  scale_fill_manual(values = c("#0000a7","#c1272d"))+
  geom_density_ridges(aes(y=outcome),
                      jittered_points = TRUE,
                      position = "raincloud",
                      alpha = 0.7, scale = 0.5)+
  labs(x="Load (% of 1RM) Per Set",y="Density")

primary.densitites.frequency<-ggplot(data = primary.data,aes(x=frequency.per.muscle,color=outcome,fill=outcome))+THEME+
  theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank())+
  scale_x_continuous(breaks = c(1,2,3,4,5,6))+
  scale_color_manual(values = c("#0000a7","#c1272d"))+
  scale_fill_manual(values = c("#0000a7","#c1272d"))+
  geom_density_ridges(aes(y=outcome),
                      jittered_points = TRUE,
                      position = "raincloud",
                      alpha = 0.7, scale = 0.5)+
  labs(x="Direct Frequency (Sessions Per Week Per Muscle Group)",y="Density")

primary.density.plot.final<-primary.densitites.volume+primary.densitites.load+primary.densitites.frequency&
  plot_annotation(title = "Visual Summaries of Training Interventions Included in Meta-regression Models")&
  theme(plot.title = element_text(size = 15,hjust = 0.5))&THEME

ggsave(primary.density.plot.final,device = "pdf",filename="primary.density.plot.final.pdf",width = 12,height = 6,
       path = "/Users/zacrobinson/Documents/For Publication PTF Meta/Plots")

#Study Characteristics
primary.chars<-c('train.status','nutrition.controlled','formal.cardio.intervention','progression','failure.definition',
                 'alternative.set.structure','set.rep.equated','upper.lower.other','train.exercise','within.between.design','train.intent')

for (x in primary.chars) {
  list<-primary.data%>%
    group_by(outcome,.data[[x]])%>%
    count(.data[[x]])%>%
    spread(outcome,n)%>%
    mutate(Variable=x)%>%
    rename(Level=.data[[x]])%>%
    relocate(Variable)
  
  assign((paste0("primary.char.table.",x)),list)
}

primary.char.table.combined<-rbind(primary.char.table.train.status,
                                   primary.char.table.nutrition.controlled,
                                   primary.char.table.formal.cardio.intervention,
                                   primary.char.table.progression,
                                   primary.char.table.failure.definition,
                                   primary.char.table.alternative.set.structure,
                                   primary.char.table.set.rep.equated,
                                   primary.char.table.upper.lower.other,
                                   primary.char.table.train.exercise,
                                   primary.char.table.within.between.design,
                                   primary.char.table.train.intent)

primary.char.table.combined$Variable=c(rep("Training Status",2),
                                       rep("Nutrition Control",2),
                                       rep("Concurrent Training Intervention",2),
                                       rep("Employed Progressive Overload",2),
                                       rep("Failure Definition Provided",2),
                                       rep("Alternative Set Structures",2),
                                       rep("Method of Equating Volume",3),
                                       rep("Upper or Lower Body Outcome",2),
                                       rep("Training Exercise Selection",3),
                                       rep("Within or Between Participant Design",2),
                                       rep("Utilized Maximal Intended Concentric Velocity",2))

primary.char.table.combined$Level=c("Trained","Untrained",
                                    "No","Yes",
                                    "No","Yes",
                                    "No","Yes",
                                    "No","Yes",
                                    "No","Yes",
                                    "Set","Repetition","Both",
                                    "Lower","Upper","Both",
                                    "Multi-joint","Single-joint",
                                    "Between","Within",
                                    "No","Yes")

primary.char.table.combined.final<-primary.char.table.combined%>%
  gt(groupname_col =  "Variable")%>%
  tab_header("Characteristics of Effects Included in Meta-regression Models")%>%
  opt_stylize(style = 1,color = "gray")%>%
  tab_options(row_group.as_column = TRUE)%>%
  tab_style(style = cell_text(style = "italic"),locations = cells_row_groups())%>%
  tab_style(style = cell_text(align = "center"),locations = cells_column_labels())%>%
  tab_style(style = cell_text(weight  = "bold"),locations = cells_title())%>%
  tab_style(style = cell_text(align = "center"),locations = cells_body())

gtsave(primary.char.table.combined.final,filename = "primary.char.table.combined.final.png",
       path = "/Users/zacrobinson/Documents/For Publication PTF Meta/Tables")

# FIT AND COMPARE MODELS --------------------------------------------------------------

### strength
#### linear
primary.str.g.linear.model.ri=rma.mv(yi,vi,data = primary.data.str,
                                  mods = ~ avg.rir+load.set+set.rep.equated+weeks+train.status,
                                  random = list(~1|study,~1|group,~1|obs),
                                  method = "REML", test = "t",dfs = "contain")

primary.str.g.linear.model.rs=rma.mv(yi,vi,data = primary.data.str,
                                  mods = ~ avg.rir+load.set+set.rep.equated+weeks+train.status,
                                  random = list(~avg.rir|study,~1|group,~1|obs),
                                  method = "REML", test = "t",dfs = "contain",struct = "GEN")

primary.str.g.linear.model.rs2=rma.mv(yi,vi,data = primary.data.str,
                                  mods = ~ avg.rir+load.set+set.rep.equated+weeks+train.status,
                                  random = list(~avg.rir|study,~avg.rir|group,~1|obs),
                                  method = "REML", test = "t",dfs = "contain",struct = c("GEN","GEN"))

primary.str.rom.linear.model.ri=rma.mv(yi.rom,vi.rom,data = primary.data.str,
                                  mods = ~ avg.rir+load.set+set.rep.equated+weeks+train.status,
                                  random = list(~1|study,~1|group,~1|obs),
                                  method = "REML", test = "t",dfs = "contain")

primary.str.rom.linear.model.rs=rma.mv(yi.rom,vi.rom,data = primary.data.str,
                                    mods = ~ avg.rir+load.set+set.rep.equated+weeks+train.status,
                                    random = list(~avg.rir|study,~1|group,~1|obs),
                                    method = "REML", test = "t",dfs = "contain",struct = "GEN")

primary.str.rom.linear.model.rs2=rma.mv(yi.rom,vi.rom,data = primary.data.str,
                                    mods = ~ avg.rir+load.set+set.rep.equated+weeks+train.status,
                                    random = list(~avg.rir|study,~avg.rir|group,~1|obs),
                                    method = "REML", test = "t",dfs = "contain",struct = c("GEN","GEN"))

#### linear-log
primary.str.g.log.model.ri=rma.mv(yi,vi,data = primary.data.str,
                               mods = ~ log1p(avg.rir)+load.set+set.rep.equated+weeks+train.status,
                               random = list(~1|study,~1|group,~1|obs),
                               method = "REML", test = "t",dfs = "contain")

primary.str.g.log.model.rs=rma.mv(yi,vi,data = primary.data.str,
                               mods = ~ log1p(avg.rir)+load.set+set.rep.equated+weeks+train.status,
                               random = list(~avg.rir|study,~1|group,~1|obs),
                               method = "REML", test = "t",dfs = "contain",struct = "GEN")

primary.str.g.log.model.rs2=rma.mv(yi,vi,data = primary.data.str,
                                  mods = ~ log1p(avg.rir)+load.set+set.rep.equated+weeks+train.status,
                                  random = list(~avg.rir|study,~avg.rir|group,~1|obs),
                                  method = "REML", test = "t",dfs = "contain",struct = c("GEN","GEN"))

primary.str.rom.log.model.ri=rma.mv(yi.rom,vi.rom,data = primary.data.str,
                                    mods = ~ log1p(avg.rir)+load.set+set.rep.equated+weeks+train.status,
                                    random = list(~1|study,~1|group,~1|obs),
                                    method = "REML", test = "t",dfs = "contain")

primary.str.rom.log.model.rs=rma.mv(yi.rom,vi.rom,data = primary.data.str,
                                 mods = ~ log1p(avg.rir)+load.set+set.rep.equated+weeks+train.status,
                                 random = list(~avg.rir|study,~1|group,~1|obs),
                                 method = "REML", test = "t",dfs = "contain",struct = "GEN")

primary.str.rom.log.model.rs2=rma.mv(yi.rom,vi.rom,data = primary.data.str,
                                 mods = ~ log1p(avg.rir)+load.set+set.rep.equated+weeks+train.status,
                                 random = list(~avg.rir|study,~avg.rir|group,~1|obs),
                                 method = "REML", test = "t",dfs = "contain",struct = c("GEN","GEN"))

#### quadratic
primary.str.g.quadratic.model.ri=rma.mv(yi,vi,data = primary.data.str,
                                  mods = ~ poly(avg.rir,2,raw = TRUE)+load.set+set.rep.equated+weeks+train.status,
                                  random = list(~1|study,~1|group,~1|obs),
                                  method = "REML", test = "t",dfs = "contain")

primary.str.g.quadratic.model.rs=rma.mv(yi,vi,data = primary.data.str,
                                     mods = ~ poly(avg.rir,2,raw = TRUE)+load.set+set.rep.equated+weeks+train.status,
                                     random = list(~avg.rir|study,~1|group,~1|obs),
                                     method = "REML", test = "t",dfs = "contain",struct = "GEN")

primary.str.g.quadratic.model.rs2=rma.mv(yi,vi,data = primary.data.str,
                                     mods = ~ poly(avg.rir,2,raw = TRUE)+load.set+set.rep.equated+weeks+train.status,
                                     random = list(~avg.rir|study,~avg.rir|group,~1|obs),
                                     method = "REML", test = "t",dfs = "contain",struct = c("GEN","GEN"))

primary.str.rom.quadratic.model.ri=rma.mv(yi.rom,vi.rom,data = primary.data.str,
                                    mods = ~ poly(avg.rir,2,raw = TRUE)+load.set+set.rep.equated+weeks+train.status,
                                    random = list(~1|study,~1|group,~1|obs),
                                    method = "REML", test = "t",dfs = "contain")

primary.str.rom.quadratic.model.rs=rma.mv(yi.rom,vi.rom,data = primary.data.str,
                                       mods = ~ poly(avg.rir,2,raw = TRUE)+load.set+set.rep.equated+weeks+train.status,
                                       random = list(~avg.rir|study,~1|group,~1|obs),
                                       method = "REML", test = "t",dfs = "contain",struct = "GEN")

primary.str.rom.quadratic.model.rs2=rma.mv(yi.rom,vi.rom,data = primary.data.str,
                                       mods = ~ poly(avg.rir,2,raw = TRUE)+load.set+set.rep.equated+weeks+train.status,
                                       random = list(~avg.rir|study,~avg.rir|group,~1|obs),
                                       method = "REML", test = "t",dfs = "contain",struct = c("GEN","GEN"))

#### cubic
primary.str.g.cubic.model.ri=rma.mv(yi,vi,data = primary.data.str,
                                     mods = ~ poly(avg.rir,3,raw = TRUE)+load.set+set.rep.equated+weeks+train.status,
                                     random = list(~1|study,~1|group,~1|obs),
                                 method = "REML", test = "t",dfs = "contain")

primary.str.g.cubic.model.rs=rma.mv(yi,vi,data = primary.data.str,
                                 mods = ~ poly(avg.rir,3,raw = TRUE)+load.set+set.rep.equated+weeks+train.status,
                                 random = list(~avg.rir|study,~1|group,~1|obs),
                                 method = "REML", test = "t",dfs = "contain",struct = "GEN")

primary.str.g.cubic.model.rs2=rma.mv(yi,vi,data = primary.data.str,
                                 mods = ~ poly(avg.rir,3,raw = TRUE)+load.set+set.rep.equated+weeks+train.status,
                                 random = list(~avg.rir|study,~avg.rir|group,~1|obs),
                                 method = "REML", test = "t",dfs = "contain",struct = c("GEN","GEN"))

primary.str.rom.cubic.model.ri=rma.mv(yi.rom,vi.rom,data = primary.data.str,
                                       mods = ~ poly(avg.rir,3,raw = TRUE)+load.set+set.rep.equated+weeks+train.status,
                                       random = list(~1|study,~1|group,~1|obs),
                                   method = "REML", test = "t",dfs = "contain")

primary.str.rom.cubic.model.rs=rma.mv(yi.rom,vi.rom,data = primary.data.str,
                                   mods = ~ poly(avg.rir,3,raw = TRUE)+load.set+set.rep.equated+weeks+train.status,
                                   random = list(~avg.rir|study,~1|group,~1|obs),
                                   method = "REML", test = "t",dfs = "contain",struct = "GEN")

primary.str.rom.cubic.model.rs2=rma.mv(yi.rom,vi.rom,data = primary.data.str,
                                   mods = ~ poly(avg.rir,3,raw = TRUE)+load.set+set.rep.equated+weeks+train.status,
                                   random = list(~avg.rir|study,~avg.rir|group,~1|obs),
                                   method = "REML", test = "t",dfs = "contain",struct = c("GEN","GEN"))

#### linear-spline
primary.str.g.spline.model.ri=rma.mv(yi,vi,data = primary.data.str,
                                  mods = ~ bs(spline.rir,degree = 1,knots = 0)+load.set+set.rep.equated+weeks+train.status,
                                  random = list(~1|study,~1|group,~1|obs),
                                  method = "REML", test = "t",dfs = "contain")

primary.str.g.spline.model.rs=rma.mv(yi,vi,data = primary.data.str,
                                  mods = ~ bs(spline.rir,degree = 1,knots = 0)+load.set+set.rep.equated+weeks+train.status,
                                  random = list(~spline.rir|study,~1|group,~1|obs),
                                  method = "REML", test = "t",dfs = "contain",struct = "GEN")

primary.str.g.spline.model.rs2=rma.mv(yi,vi,data = primary.data.str,
                                  mods = ~ bs(spline.rir,degree = 1,knots = 0)+load.set+set.rep.equated+weeks+train.status,
                                  random = list(~spline.rir|study,~spline.rir|group,~1|obs),
                                  method = "REML", test = "t",dfs = "contain",struct = c("GEN","GEN"))

primary.str.rom.spline.model.ri=rma.mv(yi.rom,vi.rom,data = primary.data.str,
                                    mods = ~ bs(spline.rir,degree = 1,knots = 0)+load.set+set.rep.equated+weeks+train.status,
                                    random = list(~1|study,~1|group,~1|obs),
                                    method = "REML", test = "t",dfs = "contain")

primary.str.rom.spline.model.rs=rma.mv(yi.rom,vi.rom,data = primary.data.str,
                                    mods = ~ bs(spline.rir,degree = 1,knots = 0)+load.set+set.rep.equated+weeks+train.status,
                                    random = list(~spline.rir|study,~1|group,~1|obs),
                                    method = "REML", test = "t",dfs = "contain",struct = "GEN")

primary.str.rom.spline.model.rs2=rma.mv(yi.rom,vi.rom,data = primary.data.str,
                                    mods = ~ bs(spline.rir,degree = 1,knots = 0)+load.set+set.rep.equated+weeks+train.status,
                                    random = list(~spline.rir|study,~spline.rir|group,~1|obs),
                                    method = "REML", test = "t",dfs = "contain",struct = c("GEN","GEN"))

#### cubic-spline
primary.str.g.rcs.model.ri=rma.mv(yi,vi,data = primary.data.str,
                                  mods = ~ ns(avg.rir,knots = primary.str.k,Boundary.knots=primary.str.b)+load.set+set.rep.equated+weeks+train.status,
                                  random = list(~1|study,~1|group,~1|obs),
                               method = "REML", test = "t",dfs = "contain")

primary.str.g.rcs.model.rs=rma.mv(yi,vi,data = primary.data.str,
                               mods = ~ ns(avg.rir,knots = primary.str.k,Boundary.knots = primary.str.b)+load.set+set.rep.equated+weeks+train.status,
                               random = list(~avg.rir|study,~1|group,~1|obs),
                               method = "REML", test = "t",dfs = "contain",struct = "GEN")

primary.str.g.rcs.model.rs2=rma.mv(yi,vi,data = primary.data.str,
                               mods = ~ ns(avg.rir,knots = primary.str.k,Boundary.knots = primary.str.b)+load.set+set.rep.equated+weeks+train.status,
                               random = list(~avg.rir|study,~avg.rir|group,~1|obs),
                               method = "REML", test = "t",dfs = "contain",struct = c("GEN","GEN"))

primary.str.rom.rcs.model.ri=rma.mv(yi.rom,vi.rom,data = primary.data.str,
                                    mods = ~ns(avg.rir,knots = primary.str.k,Boundary.knots = primary.str.b)+load.set+set.rep.equated+weeks+train.status,
                                    random = list(~1|study,~1|group,~1|obs),
                                 method = "REML", test = "t",dfs = "contain")

primary.str.rom.rcs.model.rs=rma.mv(yi.rom,vi.rom,data = primary.data.str,
                                 mods = ~ns(avg.rir,knots = primary.str.k,Boundary.knots = primary.str.b)+load.set+set.rep.equated+weeks+train.status,
                                 random = list(~avg.rir|study,~1|group,~1|obs),
                                 method = "REML", test = "t",dfs = "contain",struct = "GEN")

primary.str.rom.rcs.model.rs2=rma.mv(yi.rom,vi.rom,data = primary.data.str,
                                 mods = ~ns(avg.rir,knots = primary.str.k,Boundary.knots = primary.str.b)+load.set+set.rep.equated+weeks+train.status,
                                 random = list(~avg.rir|study,~avg.rir|group,~1|obs),
                                 method = "REML", test = "t",dfs = "contain",struct = c("GEN","GEN"))

#### compare g models
primary.model.compare.str.g=bf_models(primary.str.g.linear.model.ri,
          primary.str.g.linear.model.rs,
          primary.str.g.linear.model.rs2,
          primary.str.g.log.model.ri,
          primary.str.g.log.model.rs,
          primary.str.g.log.model.rs2,
          primary.str.g.quadratic.model.ri,
          primary.str.g.quadratic.model.rs,
          primary.str.g.quadratic.model.rs2,
          primary.str.g.cubic.model.ri,
          primary.str.g.cubic.model.rs,
          primary.str.g.cubic.model.rs2,
          primary.str.g.spline.model.ri,
          primary.str.g.spline.model.rs,
          primary.str.g.spline.model.rs2,
          primary.str.g.rcs.model.ri,
          primary.str.g.rcs.model.rs,
          primary.str.g.rcs.model.rs2
          )%>%as.matrix()%>%as.data.frame()%>%
  rename(linear.ri=1,
         linear.rs=2,
         linear.rs2=3,
         log.ri=4,
         log.rs=5,
         log.rs2=6,
         quadratic.ri=7,
         quadratic.rs=8,
         quadratic.rs2=9,
         cubic.ri=10,
         cubic.rs=11,
         cubic.rs2=12,
         spline.ri=13,
         spline.rs=14,
         spline.rs2=15,
         rcs.ri=16,
         rcs.rs=17,
         rcs.rs2=18)%>%
  mutate(model=c(
    "linear.ri",
    "linear.rs",
    "linear.rs2",
    "log.ri",
    "log.rs",
    "log.rs2",
    "quadratic.ri",
    "quadratic.rs",
    "quadratic.rs2",
    "cubic.ri",
    "cubic.rs",
    "cubic.rs2",
    "spline.ri",
    "spline.rs",
    "spline.rs2",
    "rcs.ri",
    "rcs.rs",
    "rcs.rs2"))%>%
  remove_rownames()%>%
  melt()%>%
  mutate(value=round(2*value,2))%>%
ggplot(aes(x=fct_inorder(variable),
           y=fct_relevel(model,
                         "rcs.rs2",
                         "rcs.rs",
                         "rcs.ri",
                         "spline.rs2",
                         "spline.rs",
                         "spline.ri",
                         "cubic.rs2",
                         "cubic.rs",
                         "cubic.ri",
                         "quadratic.rs2",
                         "quadratic.rs",
                         "quadratic.ri",
                         "log.rs2",
                         "log.rs",
                         "log.ri",
                         "linear.rs2",
                         "linear.rs",
                         "linear.ri")))+THEME+theme(legend.position = "right",
                                              plot.subtitle = element_text(face = c("italic")))+
  labs(x="Numerator",
       y="Denominator",
       title = "Maximal Strength Model Comparision Matrix (Standardized Mean Change)",
       subtitle = "2 x Log(Bayes Factor) Positive Values Favor Numerator",
       fill="2 x Log(BF)",
       caption = "Kass and Raferty (1995) Scale: -Inf to 0 = Negative; 0 to 2 = Weak; 2 to 6 = Positive; 6 to 10 = Strong; 10 to +Inf = Very Strong")+
  geom_tile(aes(fill=value))+
  geom_text(aes(label=value),color="black")+
  scale_fill_gradient2(low = "gray50",
                      high = "#c1272d",
                      mid="white",
                      midpoint = 0)+
  scale_x_discrete(labels=c("Linear\n(Intercept)",
                            "Linear\n(Slope)",
                            "Linear\n(Slopes)",
                            "Linear-log\n(Intercept)",
                            "Linear-log\n(Slope)",
                            "Linear-log\n(Slopes)",
                            "Quadratic\n(Intercept)",
                            "Quadratic\n(Slope)",
                            "Quadratic\n(Slopes)",
                            "Cubic\n(Intercept)",
                            "Cubic\n(Slope)",
                            "Cubic\n(Slopes)",
                            "Linear\nSpline\n(Intercept)",
                            "Linear\nSpline\n(Slope)",
                            "Linear\nSpline\n(Slopes)",
                            "Restricted\nCubic Spline\n(Intercept)",
                            "Restricted\nCubic Spline\n(Slope)",
                            "Restricted\nCubic Spline\n(Slopes)"),
                   position = "top")+
  scale_y_discrete(labels=rev(c("Linear\n(Intercept)",
                            "Linear\n(Slope)",
                            "Linear\n(Slopes)",
                            "Linear-log\n(Intercept)",
                            "Linear-log\n(Slope)",
                            "Linear-log\n(Slopes)",
                            "Quadratic\n(Intercept)",
                            "Quadratic\n(Slope)",
                            "Quadratic\n(Slopes)",
                            "Cubic\n(Intercept)",
                            "Cubic\n(Slope)",
                            "Cubic\n(Slopes)",
                            "Linear\nSpline\n(Intercept)",
                            "Linear\nSpline\n(Slope)",
                            "Linear\nSpline\n(Slopes)",
                            "Restricted\nCubic Spline\n(Intercept)",
                            "Restricted\nCubic Spline\n(Slope)", 
                            "Restricted\nCubic Spline\n(Slopes)"))) # log.ri best fit

ggsave(primary.model.compare.str.g,device = "pdf",filename="primary.model.compare.str.g.pdf",width = 14,height = 9,
       path = "/Users/zacrobinson/Documents/For Publication PTF Meta/Plots")

#### refit g model with robust variance estimates
primary.str.g.best.model.prep=rma.mv(yi,vi,data = primary.data.str,
                                  mods = ~ log1p(avg.rir)+load.set+set.rep.equated+weeks+train.status,
                                  random = list(~1|study,~1|group,~1|obs),
                                  method = "REML", test = "t",dfs = "contain")

primary.str.g.best.model=robust.rma.mv(primary.str.g.best.model.prep,
                                       cluster = primary.data.str$study)

#### I2 and R2 for g model
primary.str.g.best.model.r2=r2_ml(primary.str.g.best.model)*100
save(primary.str.g.best.model.r2,file = "primary.str.g.best.model.r2.RData")
primary.str.g.best.model.i2=i2_ml(primary.str.g.best.model)
save(primary.str.g.best.model.i2,file = "primary.str.g.best.model.i2.RData")

#### g model output
primary.str.g.best.output<-subset(as.data.frame(coef(summary(primary.str.g.best.model))),select = c(-pval,-tval))
primary.str.g.best.output.prep1<-qt(.975,primary.str.g.best.output$df)
primary.str.g.best.output.prep2<-sqrt((primary.str.g.best.output$se^2)+sum(primary.str.g.best.model$sigma2))
primary.str.g.best.output$lower.PI<-primary.str.g.best.output$estimate-(primary.str.g.best.output.prep1*primary.str.g.best.output.prep2)
primary.str.g.best.output$upper.PI<-primary.str.g.best.output$estimate+(primary.str.g.best.output.prep1*primary.str.g.best.output.prep2)
primary.str.g.best.output$Parameter<-c("Intercept",
                                       "Log<sub>+1</sub>(Estimated RIR)",
                                       "Load",
                                       "Volume Equated [Rep]",
                                       "Volume Equated [Both]",
                                       "Weeks",
                                       "Training Status [Untrained]")

primary.str.g.best.output.table<-relocate(primary.str.g.best.output,Parameter)%>%
  gt()%>%
  cols_label(ci.lb="Lower Bound",
             ci.ub="Upper Bound",
             lower.PI="Lower Bound",
             upper.PI="Upper Bound",
             se="Standard Error",
             estimate="Estimate")%>%
  fmt_markdown(columns = 1)%>%
  tab_header(title = "Coefficients of Linear-log Multilevel Model for Maximal Strength (Standardized Mean Change)")%>%
  tab_spanner(columns = 5:6,label = "95% Confidence Interval")%>%
  tab_spanner(columns = 7:8,label = "95% Predicition Interval")%>%
  opt_stylize(style = 1,color = "gray")%>%
  tab_options(column_labels.background.color ="gray97",column_labels.font.weight = "bold")%>%
  tab_style(style = cell_text(style = "italic"),locations = cells_row_groups())%>%
  tab_style(style = cell_text(align = "center"),locations = cells_column_labels())%>%
  tab_style(style = cell_text(weight  = "bold"),locations = cells_title())%>%
  tab_style(style = cell_text(align = "center"),locations = cells_body(columns=-1))%>%
  cols_width(-df~px(133))%>%
  cols_width(Parameter~px(200))

gtsave(primary.str.g.best.output.table,filename = "primary.str.g.best.output.table.png",
       path = "/Users/zacrobinson/Documents/For Publication PTF Meta/Models")

#### compare rom models
primary.model.compare.str.rom=bf_models(primary.str.rom.linear.model.ri,
                                      primary.str.rom.linear.model.rs,
                                      primary.str.rom.linear.model.rs2,
                                      primary.str.rom.log.model.ri,
                                      primary.str.rom.log.model.rs,
                                      primary.str.rom.log.model.rs2,
                                      primary.str.rom.quadratic.model.ri,
                                      primary.str.rom.quadratic.model.rs,
                                      primary.str.rom.quadratic.model.rs2,
                                      primary.str.rom.cubic.model.ri,
                                      primary.str.rom.cubic.model.rs,
                                      primary.str.rom.cubic.model.rs2,
                                      primary.str.rom.spline.model.ri,
                                      primary.str.rom.spline.model.rs,
                                      primary.str.rom.spline.model.rs2,
                                      primary.str.rom.rcs.model.ri,
                                      primary.str.rom.rcs.model.rs,
                                      primary.str.rom.rcs.model.rs2
)%>%as.matrix()%>%as.data.frame()%>%
  rename(linear.ri=1,
         linear.rs=2,
         linear.rs2=3,
         log.ri=4,
         log.rs=5,
         log.rs2=6,
         quadratic.ri=7,
         quadratic.rs=8,
         quadratic.rs2=9,
         cubic.ri=10,
         cubic.rs=11,
         cubic.rs2=12,
         spline.ri=13,
         spline.rs=14,
         spline.rs2=15,
         rcs.ri=16,
         rcs.rs=17,
         rcs.rs2=18)%>%
  mutate(model=c(
    "linear.ri",
    "linear.rs",
    "linear.rs2",
    "log.ri",
    "log.rs",
    "log.rs2",
    "quadratic.ri",
    "quadratic.rs",
    "quadratic.rs2",
    "cubic.ri",
    "cubic.rs",
    "cubic.rs2",
    "spline.ri",
    "spline.rs",
    "spline.rs2",
    "rcs.ri",
    "rcs.rs",
    "rcs.rs2"))%>%
  remove_rownames()%>%
  melt()%>%
  mutate(value=round(2*value,2))%>% 
  ggplot(aes(x=fct_inorder(variable),
             y=fct_relevel(model,
                           "rcs.rs2",
                           "rcs.rs",
                           "rcs.ri",
                           "spline.rs2",
                           "spline.rs",
                           "spline.ri",
                           "cubic.rs2",
                           "cubic.rs",
                           "cubic.ri",
                           "quadratic.rs2",
                           "quadratic.rs",
                           "quadratic.ri",
                           "log.rs2",
                           "log.rs",
                           "log.ri",
                           "linear.rs2",
                           "linear.rs",
                           "linear.ri")))+THEME+theme(legend.position = "right",
                                                      plot.subtitle = element_text(face = c("italic")))+
  labs(x="Numerator",
       y="Denominator",
       title = "Maximal Strength Model Comparision Matrix (Response Ratio)",
       subtitle = "2 x Log(Bayes Factor) Positive Values Favor Numerator",
       fill="2 x Log(BF)",
       caption = "Kass and Raferty (1995) Scale: -Inf to 0 = Negative; 0 to 2 = Weak; 2 to 6 = Positive; 6 to 10 = Strong; 10 to +Inf = Very Strong")+
  geom_tile(aes(fill=value))+
  geom_text(aes(label=value),color="black")+
  scale_fill_gradient2(low = "gray50",
                       high = "#c1272d",
                       mid="white",
                       midpoint = 0)+
  scale_x_discrete(labels=c("Linear\n(Intercept)",
                            "Linear\n(Slope)",
                            "Linear\n(Slopes)",
                            "Linear-log\n(Intercept)",
                            "Linear-log\n(Slope)",
                            "Linear-log\n(Slopes)",
                            "Quadratic\n(Intercept)",
                            "Quadratic\n(Slope)",
                            "Quadratic\n(Slopes)",
                            "Cubic\n(Intercept)",
                            "Cubic\n(Slope)",
                            "Cubic\n(Slopes)",
                            "Linear\nSpline\n(Intercept)",
                            "Linear\nSpline\n(Slope)",
                            "Linear\nSpline\n(Slopes)",
                            "Restricted\nCubic Spline\n(Intercept)",
                            "Restricted\nCubic Spline\n(Slope)",
                            "Restricted\nCubic Spline\n(Slopes)"))+
  scale_y_discrete(labels=rev(c("Linear\n(Intercept)",
                                "Linear\n(Slope)",
                                "Linear\n(Slopes)",
                                "Linear-log\n(Intercept)",
                                "Linear-log\n(Slope)",
                                "Linear-log\n(Slopes)",
                                "Quadratic\n(Intercept)",
                                "Quadratic\n(Slope)",
                                "Quadratic\n(Slopes)",
                                "Cubic\n(Intercept)",
                                "Cubic\n(Slope)",
                                "Cubic\n(Slopes)",
                                "Linear\nSpline\n(Intercept)",
                                "Linear\nSpline\n(Slope)",
                                "Linear\nSpline\n(Slopes)",
                                "Restricted\nCubic Spline\n(Intercept)",
                                "Restricted\nCubic Spline\n(Slope)", 
                                "Restricted\nCubic Spline\n(Slopes)"))) # linear.ri best fit

ggsave(primary.model.compare.str.rom,device = "pdf",filename="primary.model.compare.str.rom.pdf",width = 14,height = 9,
       path = "/Users/zacrobinson/Documents/For Publication PTF Meta/Plots")

#### refit rom model with robust variance estimates
primary.str.rom.best.model.prep=rma.mv(yi.rom,vi.rom,data = primary.data.str,
                                mods = ~avg.rir+load.set+set.rep.equated+weeks+train.status,
                                random = list(~1|study,~1|group,~1|obs),
                                method = "REML", test = "t",dfs = "contain")

primary.str.rom.best.model=robust.rma.mv(primary.str.rom.best.model.prep,
                                       cluster = primary.data.str$study)

#### I2 and R2 for rom model
primary.str.rom.best.model.r2=r2_ml(primary.str.rom.best.model)*100
save(primary.str.rom.best.model.r2,file = "primary.str.rom.best.model.r2.RData")
primary.str.rom.best.model.i2=i2_ml(primary.str.rom.best.model)
save(primary.str.rom.best.model.i2,file = "primary.str.rom.best.model.i2.RData")

#### rom model output
primary.str.rom.best.output<-subset(as.data.frame(coef(summary(primary.str.rom.best.model))),select = c(-pval,-tval))
primary.str.rom.best.output.prep1<-qt(.975,primary.str.rom.best.output$df)
primary.str.rom.best.output.prep2<-sqrt((primary.str.rom.best.output$se^2)+sum(primary.str.rom.best.model$sigma2))
primary.str.rom.best.output$lower.PI<-primary.str.rom.best.output$estimate-(primary.str.rom.best.output.prep1*primary.str.rom.best.output.prep2)
primary.str.rom.best.output$upper.PI<-primary.str.rom.best.output$estimate+(primary.str.rom.best.output.prep1*primary.str.rom.best.output.prep2)
primary.str.rom.best.output$Parameter<-c("Intercept",
                                         "Estimated RIR",
                                         "Load",
                                         "Volume Equated [Rep]",
                                         "Volume Equated [Both]",
                                         "Weeks",
                                         "Training Status [Untrained]")

primary.str.rom.best.output.table<-relocate(primary.str.rom.best.output,Parameter)%>%
  gt()%>%
  cols_label(ci.lb="Lower Bound",
             ci.ub="Upper Bound",
             lower.PI="Lower Bound",
             upper.PI="Upper Bound",
             se="Standard Error",
             estimate="Estimate")%>%
  fmt_markdown(columns = 1)%>%
  tab_header(title = "Coefficients of Linear Multilevel Model for Maximal Strength (Response Ratio)")%>%
  tab_spanner(columns = 5:6,label = "95% Confidence Interval")%>%
  tab_spanner(columns = 7:8,label = "95% Predicition Interval")%>%
  opt_stylize(style = 1,color = "gray")%>%
  tab_options(column_labels.background.color ="gray97",column_labels.font.weight = "bold")%>%
  tab_style(style = cell_text(style = "italic"),locations = cells_row_groups())%>%
  tab_style(style = cell_text(align = "center"),locations = cells_column_labels())%>%
  tab_style(style = cell_text(weight  = "bold"),locations = cells_title())%>%
  tab_style(style = cell_text(align = "center"),locations = cells_body(columns=-1))%>%
  cols_width(-df~px(133))%>%
  cols_width(Parameter~px(200))

gtsave(primary.str.rom.best.output.table,filename = "primary.str.rom.best.output.table.png",
       path = "/Users/zacrobinson/Documents/For Publication PTF Meta/Models")

### hypertrophy
#### linear
primary.hyp.g.linear.model.ri=rma.mv(yi,vi,data = primary.data.hyp,
                                     mods = ~ avg.rir+load.set+set.rep.equated+weeks+train.status,
                                     random = list(~1|study,~1|group,~1|obs),
                                     method = "REML", test = "t",dfs = "contain")

primary.hyp.g.linear.model.rs=rma.mv(yi,vi,data = primary.data.hyp,
                                     mods = ~ avg.rir+load.set+set.rep.equated+weeks+train.status,
                                     random = list(~avg.rir|study,~1|group,~1|obs),
                                     method = "REML", test = "t",dfs = "contain",struct = "GEN")

primary.hyp.g.linear.model.rs2=rma.mv(yi,vi,data = primary.data.hyp,
                                      mods = ~ avg.rir+load.set+set.rep.equated+weeks+train.status,
                                      random = list(~avg.rir|study,~avg.rir|group,~1|obs),
                                      method = "REML", test = "t",dfs = "contain",struct = c("GEN","GEN"))

primary.hyp.rom.linear.model.ri=rma.mv(yi.rom,vi.rom,data = primary.data.hyp,
                                       mods = ~ avg.rir+load.set+set.rep.equated+weeks+train.status,
                                       random = list(~1|study,~1|group,~1|obs),
                                       method = "REML", test = "t",dfs = "contain")

primary.hyp.rom.linear.model.rs=rma.mv(yi.rom,vi.rom,data = primary.data.hyp,
                                       mods = ~ avg.rir+load.set+set.rep.equated+weeks+train.status,
                                       random = list(~avg.rir|study,~1|group,~1|obs),
                                       method = "REML", test = "t",dfs = "contain",struct = "GEN")

primary.hyp.rom.linear.model.rs2=rma.mv(yi.rom,vi.rom,data = primary.data.hyp,
                                        mods = ~ avg.rir+load.set+set.rep.equated+weeks+train.status,
                                        random = list(~avg.rir|study,~avg.rir|group,~1|obs),
                                        method = "REML", test = "t",dfs = "contain",struct = c("GEN","GEN"))

#### linear-log
primary.hyp.g.log.model.ri=rma.mv(yi,vi,data = primary.data.hyp,
                                  mods = ~ log1p(avg.rir)+load.set+set.rep.equated+weeks+train.status,
                                  random = list(~1|study,~1|group,~1|obs),
                                  method = "REML", test = "t",dfs = "contain")

primary.hyp.g.log.model.rs=rma.mv(yi,vi,data = primary.data.hyp,
                                  mods = ~ log1p(avg.rir)+load.set+set.rep.equated+weeks+train.status,
                                  random = list(~avg.rir|study,~1|group,~1|obs),
                                  method = "REML", test = "t",dfs = "contain",struct = "GEN")

primary.hyp.g.log.model.rs2=rma.mv(yi,vi,data = primary.data.hyp,
                                   mods = ~ log1p(avg.rir)+load.set+set.rep.equated+weeks+train.status,
                                   random = list(~avg.rir|study,~avg.rir|group,~1|obs),
                                   method = "REML", test = "t",dfs = "contain",struct = c("GEN","GEN"))

primary.hyp.rom.log.model.ri=rma.mv(yi.rom,vi.rom,data = primary.data.hyp,
                                    mods = ~ log1p(avg.rir)+load.set+set.rep.equated+weeks+train.status,
                                    random = list(~1|study,~1|group,~1|obs),
                                    method = "REML", test = "t",dfs = "contain")

primary.hyp.rom.log.model.rs=rma.mv(yi.rom,vi.rom,data = primary.data.hyp,
                                    mods = ~ log1p(avg.rir)+load.set+set.rep.equated+weeks+train.status,
                                    random = list(~avg.rir|study,~1|group,~1|obs),
                                    method = "REML", test = "t",dfs = "contain",struct = "GEN")

primary.hyp.rom.log.model.rs2=rma.mv(yi.rom,vi.rom,data = primary.data.hyp,
                                     mods = ~ log1p(avg.rir)+load.set+set.rep.equated+weeks+train.status,
                                     random = list(~avg.rir|study,~avg.rir|group,~1|obs),
                                     method = "REML", test = "t",dfs = "contain",struct = c("GEN","GEN"))

#### quadratic
primary.hyp.g.quadratic.model.ri=rma.mv(yi,vi,data = primary.data.hyp,
                                        mods = ~ poly(avg.rir,2,raw = TRUE)+load.set+set.rep.equated+weeks+train.status,
                                        random = list(~1|study,~1|group,~1|obs),
                                        method = "REML", test = "t",dfs = "contain")

primary.hyp.g.quadratic.model.rs=rma.mv(yi,vi,data = primary.data.hyp,
                                        mods = ~ poly(avg.rir,2,raw = TRUE)+load.set+set.rep.equated+weeks+train.status,
                                        random = list(~avg.rir|study,~1|group,~1|obs),
                                        method = "REML", test = "t",dfs = "contain",struct = "GEN")

primary.hyp.g.quadratic.model.rs2=rma.mv(yi,vi,data = primary.data.hyp,
                                         mods = ~ poly(avg.rir,2,raw = TRUE)+load.set+set.rep.equated+weeks+train.status,
                                         random = list(~avg.rir|study,~avg.rir|group,~1|obs),
                                         method = "REML", test = "t",dfs = "contain",struct = c("GEN","GEN"))

primary.hyp.rom.quadratic.model.ri=rma.mv(yi.rom,vi.rom,data = primary.data.hyp,
                                          mods = ~ poly(avg.rir,2,raw = TRUE)+load.set+set.rep.equated+weeks+train.status,
                                          random = list(~1|study,~1|group,~1|obs),
                                          method = "REML", test = "t",dfs = "contain")

primary.hyp.rom.quadratic.model.rs=rma.mv(yi.rom,vi.rom,data = primary.data.hyp,
                                          mods = ~ poly(avg.rir,2,raw = TRUE)+load.set+set.rep.equated+weeks+train.status,
                                          random = list(~avg.rir|study,~1|group,~1|obs),
                                          method = "REML", test = "t",dfs = "contain",struct = "GEN")

primary.hyp.rom.quadratic.model.rs2=rma.mv(yi.rom,vi.rom,data = primary.data.hyp,
                                           mods = ~ poly(avg.rir,2,raw = TRUE)+load.set+set.rep.equated+weeks+train.status,
                                           random = list(~avg.rir|study,~avg.rir|group,~1|obs),
                                           method = "REML", test = "t",dfs = "contain",struct = c("GEN","GEN"))

#### cubic
primary.hyp.g.cubic.model.ri=rma.mv(yi,vi,data = primary.data.hyp,
                                    mods = ~ poly(avg.rir,3,raw = TRUE)+load.set+set.rep.equated+weeks+train.status,
                                    random = list(~1|study,~1|group,~1|obs),
                                    method = "REML", test = "t",dfs = "contain")

primary.hyp.g.cubic.model.rs=rma.mv(yi,vi,data = primary.data.hyp,
                                    mods = ~ poly(avg.rir,3,raw = TRUE)+load.set+set.rep.equated+weeks+train.status,
                                    random = list(~avg.rir|study,~1|group,~1|obs),
                                    method = "REML", test = "t",dfs = "contain",struct = "GEN")

primary.hyp.g.cubic.model.rs2=rma.mv(yi,vi,data = primary.data.hyp,
                                     mods = ~ poly(avg.rir,3,raw = TRUE)+load.set+set.rep.equated+weeks+train.status,
                                     random = list(~avg.rir|study,~avg.rir|group,~1|obs),
                                     method = "REML", test = "t",dfs = "contain",struct = c("GEN","GEN"))

primary.hyp.rom.cubic.model.ri=rma.mv(yi.rom,vi.rom,data = primary.data.hyp,
                                      mods = ~ poly(avg.rir,3,raw = TRUE)+load.set+set.rep.equated+weeks+train.status,
                                      random = list(~1|study,~1|group,~1|obs),
                                      method = "REML", test = "t",dfs = "contain")

primary.hyp.rom.cubic.model.rs=rma.mv(yi.rom,vi.rom,data = primary.data.hyp,
                                      mods = ~ poly(avg.rir,3,raw = TRUE)+load.set+set.rep.equated+weeks+train.status,
                                      random = list(~avg.rir|study,~1|group,~1|obs),
                                      method = "REML", test = "t",dfs = "contain",struct = "GEN")

primary.hyp.rom.cubic.model.rs2=rma.mv(yi.rom,vi.rom,data = primary.data.hyp,
                                       mods = ~ poly(avg.rir,3,raw = TRUE)+load.set+set.rep.equated+weeks+train.status,
                                       random = list(~avg.rir|study,~avg.rir|group,~1|obs),
                                       method = "REML", test = "t",dfs = "contain",struct = c("GEN","GEN"))

#### linear-spline
primary.hyp.g.spline.model.ri=rma.mv(yi,vi,data = primary.data.hyp,
                                     mods = ~ bs(spline.rir,degree = 1,knots = 0)+load.set+set.rep.equated+weeks+train.status,
                                     random = list(~1|study,~1|group,~1|obs),
                                     method = "REML", test = "t",dfs = "contain")

primary.hyp.g.spline.model.rs=rma.mv(yi,vi,data = primary.data.hyp,
                                     mods = ~ bs(spline.rir,degree = 1,knots = 0)+load.set+set.rep.equated+weeks+train.status,
                                     random = list(~spline.rir|study,~1|group,~1|obs),
                                     method = "REML", test = "t",dfs = "contain",struct = "GEN")

primary.hyp.g.spline.model.rs2=rma.mv(yi,vi,data = primary.data.hyp,
                                      mods = ~ bs(spline.rir,degree = 1,knots = 0)+load.set+set.rep.equated+weeks+train.status,
                                      random = list(~spline.rir|study,~spline.rir|group,~1|obs),
                                      method = "REML", test = "t",dfs = "contain",struct = c("GEN","GEN"))

primary.hyp.rom.spline.model.ri=rma.mv(yi.rom,vi.rom,data = primary.data.hyp,
                                       mods = ~ bs(spline.rir,degree = 1,knots = 0)+load.set+set.rep.equated+weeks+train.status,
                                       random = list(~1|study,~1|group,~1|obs),
                                       method = "REML", test = "t",dfs = "contain")

primary.hyp.rom.spline.model.rs=rma.mv(yi.rom,vi.rom,data = primary.data.hyp,
                                       mods = ~ bs(spline.rir,degree = 1,knots = 0)+load.set+set.rep.equated+weeks+train.status,
                                       random = list(~spline.rir|study,~1|group,~1|obs),
                                       method = "REML", test = "t",dfs = "contain",struct = "GEN")

primary.hyp.rom.spline.model.rs2=rma.mv(yi.rom,vi.rom,data = primary.data.hyp,
                                        mods = ~ bs(spline.rir,degree = 1,knots = 0)+load.set+set.rep.equated+weeks+train.status,
                                        random = list(~spline.rir|study,~spline.rir|group,~1|obs),
                                        method = "REML", test = "t",dfs = "contain",struct = c("GEN","GEN"))

#### cubic-spline
primary.hyp.g.rcs.model.ri=rma.mv(yi,vi,data = primary.data.hyp,
                                  mods = ~ ns(avg.rir,knots = primary.hyp.k,Boundary.knots = primary.hyp.b)+load.set+set.rep.equated+weeks+train.status,
                                  random = list(~1|study,~1|group,~1|obs),
                                  method = "REML", test = "t",dfs = "contain")

primary.hyp.g.rcs.model.rs=rma.mv(yi,vi,data = primary.data.hyp,
                                  mods = ~ ns(avg.rir,knots = primary.hyp.k,Boundary.knots = primary.hyp.b)+load.set+set.rep.equated+weeks+train.status,
                                  random = list(~avg.rir|study,~1|group,~1|obs),
                                  method = "REML", test = "t",dfs = "contain",struct = "GEN")

primary.hyp.g.rcs.model.rs2=rma.mv(yi,vi,data = primary.data.hyp,
                                   mods = ~ ns(avg.rir,knots = primary.hyp.k,Boundary.knots = primary.hyp.b)+load.set+set.rep.equated+weeks+train.status,
                                   random = list(~avg.rir|study,~avg.rir|group,~1|obs),
                                   method = "REML", test = "t",dfs = "contain",struct = c("GEN","GEN"))

primary.hyp.rom.rcs.model.ri=rma.mv(yi.rom,vi.rom,data = primary.data.hyp,
                                    mods = ~ns(avg.rir,knots = primary.hyp.k,Boundary.knots = primary.hyp.b)+load.set+set.rep.equated+weeks+train.status,
                                    random = list(~1|study,~1|group,~1|obs),
                                    method = "REML", test = "t",dfs = "contain")

primary.hyp.rom.rcs.model.rs=rma.mv(yi.rom,vi.rom,data = primary.data.hyp,
                                    mods = ~ns(avg.rir,knots = primary.hyp.k,Boundary.knots = primary.hyp.b)+load.set+set.rep.equated+weeks+train.status,
                                    random = list(~avg.rir|study,~1|group,~1|obs),
                                    method = "REML", test = "t",dfs = "contain",struct = "GEN")

primary.hyp.rom.rcs.model.rs2=rma.mv(yi.rom,vi.rom,data = primary.data.hyp,
                                     mods = ~ns(avg.rir,knots = primary.hyp.k,Boundary.knots = primary.hyp.b)+load.set+set.rep.equated+weeks+train.status,
                                     random = list(~avg.rir|study,~avg.rir|group,~1|obs),
                                     method = "REML", test = "t",dfs = "contain",struct = c("GEN","GEN"))

#### compare models
primary.model.compare.hyp.g=bf_models(primary.hyp.g.linear.model.ri,
                                      primary.hyp.g.linear.model.rs,
                                      primary.hyp.g.linear.model.rs2,
                                      primary.hyp.g.log.model.ri,
                                      primary.hyp.g.log.model.rs,
                                      primary.hyp.g.log.model.rs2,
                                      primary.hyp.g.quadratic.model.ri,
                                      primary.hyp.g.quadratic.model.rs,
                                      primary.hyp.g.quadratic.model.rs2,
                                      primary.hyp.g.cubic.model.ri,
                                      primary.hyp.g.cubic.model.rs,
                                      primary.hyp.g.cubic.model.rs2,
                                      primary.hyp.g.spline.model.ri,
                                      primary.hyp.g.spline.model.rs,
                                      primary.hyp.g.spline.model.rs2,
                                      primary.hyp.g.rcs.model.ri,
                                      primary.hyp.g.rcs.model.rs,
                                      primary.hyp.g.rcs.model.rs2
)%>%as.matrix()%>%as.data.frame()%>%
  rename(linear.ri=1,
         linear.rs=2,
         linear.rs2=3,
         log.ri=4,
         log.rs=5,
         log.rs2=6,
         quadratic.ri=7,
         quadratic.rs=8,
         quadratic.rs2=9,
         cubic.ri=10,
         cubic.rs=11,
         cubic.rs2=12,
         spline.ri=13,
         spline.rs=14,
         spline.rs2=15,
         rcs.ri=16,
         rcs.rs=17,
         rcs.rs2=18)%>%
  mutate(model=c(
    "linear.ri",
    "linear.rs",
    "linear.rs2",
    "log.ri",
    "log.rs",
    "log.rs2",
    "quadratic.ri",
    "quadratic.rs",
    "quadratic.rs2",
    "cubic.ri",
    "cubic.rs",
    "cubic.rs2",
    "spline.ri",
    "spline.rs",
    "spline.rs2",
    "rcs.ri",
    "rcs.rs",
    "rcs.rs2"))%>%
  remove_rownames()%>%
  melt()%>%
  mutate(value=round(2*value,2))%>% 
  ggplot(aes(x=fct_inorder(variable),
             y=fct_relevel(model,
                           "rcs.rs2",
                           "rcs.rs",
                           "rcs.ri",
                           "spline.rs2",
                           "spline.rs",
                           "spline.ri",
                           "cubic.rs2",
                           "cubic.rs",
                           "cubic.ri",
                           "quadratic.rs2",
                           "quadratic.rs",
                           "quadratic.ri",
                           "log.rs2",
                           "log.rs",
                           "log.ri",
                           "linear.rs2",
                           "linear.rs",
                           "linear.ri")))+THEME+theme(legend.position = "right",
                                                      plot.subtitle = element_text(face = c("italic")))+
  labs(x="Numerator",
       y="Denominator",
       title = "Muscle Hypertrophy Model Comparision Matrix (Standardized Mean Change)",
       subtitle = "2 x Log(Bayes Factor) Positive Values Favor Numerator",
       fill="2 x Log(BF)",
       caption = "Kass and Raferty (1995) Scale: -Inf to 0 = Negative; 0 to 2 = Weak; 2 to 6 = Positive; 6 to 10 = Strong; 10 to +Inf = Very Strong")+
  geom_tile(aes(fill=value))+
  geom_text(aes(label=value),color="black")+
  scale_fill_gradient2(low = "gray50",
                       high = "#0000a7",
                       mid="white",
                       midpoint = 0)+
  scale_x_discrete(labels=c("Linear\n(Intercept)",
                            "Linear\n(Slope)",
                            "Linear\n(Slopes)",
                            "Linear-log\n(Intercept)",
                            "Linear-log\n(Slope)",
                            "Linear-log\n(Slopes)",
                            "Quadratic\n(Intercept)",
                            "Quadratic\n(Slope)",
                            "Quadratic\n(Slopes)",
                            "Cubic\n(Intercept)",
                            "Cubic\n(Slope)",
                            "Cubic\n(Slopes)",
                            "Linear\nSpline\n(Intercept)",
                            "Linear\nSpline\n(Slope)",
                            "Linear\nSpline\n(Slopes)",
                            "Restricted\nCubic Spline\n(Intercept)",
                            "Restricted\nCubic Spline\n(Slope)",
                            "Restricted\nCubic Spline\n(Slopes)"))+
  scale_y_discrete(labels=rev(c("Linear\n(Intercept)",
                                "Linear\n(Slope)",
                                "Linear\n(Slopes)",
                                "Linear-log\n(Intercept)",
                                "Linear-log\n(Slope)",
                                "Linear-log\n(Slopes)",
                                "Quadratic\n(Intercept)",
                                "Quadratic\n(Slope)",
                                "Quadratic\n(Slopes)",
                                "Cubic\n(Intercept)",
                                "Cubic\n(Slope)",
                                "Cubic\n(Slopes)",
                                "Linear\nSpline\n(Intercept)",
                                "Linear\nSpline\n(Slope)",
                                "Linear\nSpline\n(Slopes)",
                                "Restricted\nCubic Spline\n(Intercept)",
                                "Restricted\nCubic Spline\n(Slope)", 
                                "Restricted\nCubic Spline\n(Slopes)"))) # linear.ri best fit

ggsave(primary.model.compare.hyp.g,device = "pdf",filename="primary.model.compare.hyp.g.pdf",width = 14,height = 9,
       path = "/Users/zacrobinson/Documents/For Publication PTF Meta/Plots")

#### refit g model with robust variance estimates
primary.hyp.g.best.model.prep=rma.mv(yi,vi,data = primary.data.hyp,
                                mods = ~ avg.rir+load.set+set.rep.equated+weeks+train.status,
                                random = list(~1|study,~1|group,~1|obs),
                                method = "REML", test = "t",dfs = "contain")

primary.hyp.g.best.model=robust.rma.mv(primary.hyp.g.best.model.prep,
                                       cluster = primary.data.hyp$study)
#### I2 and R2 for g model
primary.hyp.g.best.model.r2=r2_ml(primary.hyp.g.best.model)*100
save(primary.hyp.g.best.model.r2,file = "primary.hyp.g.best.model.r2.RData")
primary.hyp.g.best.model.i2=i2_ml(primary.hyp.g.best.model)
save(primary.hyp.g.best.model.i2,file = "primary.hyp.g.best.model.i2.RData")

#### g model output
primary.hyp.g.best.output<-subset(as.data.frame(coef(summary(primary.hyp.g.best.model))),select = c(-pval,-tval))
primary.hyp.g.best.output.prep1<-qt(.975,primary.hyp.g.best.output$df)
primary.hyp.g.best.output.prep2<-sqrt((primary.hyp.g.best.output$se^2)+sum(primary.hyp.g.best.model$sigma2))
primary.hyp.g.best.output$lower.PI<-primary.hyp.g.best.output$estimate-(primary.hyp.g.best.output.prep1*primary.hyp.g.best.output.prep2)
primary.hyp.g.best.output$upper.PI<-primary.hyp.g.best.output$estimate+(primary.hyp.g.best.output.prep1*primary.hyp.g.best.output.prep2)
primary.hyp.g.best.output$Parameter<-c("Intercept",
                                       "Estimated RIR",
                                       "Load",
                                       "Volume Equated [Rep]",
                                       "Volume Equated [Both]",
                                       "Weeks",
                                       "Training Status [Untrained]")

primary.hyp.g.best.output.table<-relocate(primary.hyp.g.best.output,Parameter)%>%
  gt()%>%
  cols_label(ci.lb="Lower Bound",
             ci.ub="Upper Bound",
             lower.PI="Lower Bound",
             upper.PI="Upper Bound",
             se="Standard Error",
             estimate="Estimate")%>%
  fmt_markdown(columns = 1)%>%
  tab_header(title = "Coefficients of Linear Multilevel Model for Muscle Hypertrophy (Standardized Mean Change)")%>%
  tab_spanner(columns = 5:6,label = "95% Confidence Interval")%>%
  tab_spanner(columns = 7:8,label = "95% Predicition Interval")%>%
  opt_stylize(style = 1,color = "gray")%>%
  tab_options(column_labels.background.color ="gray97",column_labels.font.weight = "bold")%>%
  tab_style(style = cell_text(style = "italic"),locations = cells_row_groups())%>%
  tab_style(style = cell_text(align = "center"),locations = cells_column_labels())%>%
  tab_style(style = cell_text(weight  = "bold"),locations = cells_title())%>%
  tab_style(style = cell_text(align = "center"),locations = cells_body(columns=-1))%>%
  cols_width(-df~px(133))%>%
  cols_width(Parameter~px(200))

gtsave(primary.hyp.g.best.output.table,filename = "primary.hyp.g.best.output.table.png",
       path = "/Users/zacrobinson/Documents/For Publication PTF Meta/Models")

#### compare rom models
primary.model.compare.hyp.rom=bf_models(primary.hyp.rom.linear.model.ri,
                                        primary.hyp.rom.linear.model.rs,
                                        primary.hyp.rom.linear.model.rs2,
                                        primary.hyp.rom.log.model.ri,
                                        primary.hyp.rom.log.model.rs,
                                        primary.hyp.rom.log.model.rs2,
                                        primary.hyp.rom.quadratic.model.ri,
                                        primary.hyp.rom.quadratic.model.rs,
                                        primary.hyp.rom.quadratic.model.rs2,
                                        primary.hyp.rom.cubic.model.ri,
                                        primary.hyp.rom.cubic.model.rs,
                                        primary.hyp.rom.cubic.model.rs2,
                                        primary.hyp.rom.spline.model.ri,
                                        primary.hyp.rom.spline.model.rs,
                                        primary.hyp.rom.spline.model.rs2,
                                        primary.hyp.rom.rcs.model.ri,
                                        primary.hyp.rom.rcs.model.rs,
                                        primary.hyp.rom.rcs.model.rs2
)%>%as.matrix()%>%as.data.frame()%>%
  rename(linear.ri=1,
         linear.rs=2,
         linear.rs2=3,
         log.ri=4,
         log.rs=5,
         log.rs2=6,
         quadratic.ri=7,
         quadratic.rs=8,
         quadratic.rs2=9,
         cubic.ri=10,
         cubic.rs=11,
         cubic.rs2=12,
         spline.ri=13,
         spline.rs=14,
         spline.rs2=15,
         rcs.ri=16,
         rcs.rs=17,
         rcs.rs2=18)%>%
  mutate(model=c(
    "linear.ri",
    "linear.rs",
    "linear.rs2",
    "log.ri",
    "log.rs",
    "log.rs2",
    "quadratic.ri",
    "quadratic.rs",
    "quadratic.rs2",
    "cubic.ri",
    "cubic.rs",
    "cubic.rs2",
    "spline.ri",
    "spline.rs",
    "spline.rs2",
    "rcs.ri",
    "rcs.rs",
    "rcs.rs2"))%>%
  remove_rownames()%>%
  melt()%>%
  mutate(value=round(2*value,2))%>% 
  ggplot(aes(x=fct_inorder(variable),
             y=fct_relevel(model,
                           "rcs.rs2",
                           "rcs.rs",
                           "rcs.ri",
                           "spline.rs2",
                           "spline.rs",
                           "spline.ri",
                           "cubic.rs2",
                           "cubic.rs",
                           "cubic.ri",
                           "quadratic.rs2",
                           "quadratic.rs",
                           "quadratic.ri",
                           "log.rs2",
                           "log.rs",
                           "log.ri",
                           "linear.rs2",
                           "linear.rs",
                           "linear.ri")))+THEME+theme(legend.position = "right",
                                                      plot.subtitle = element_text(face = c("italic")))+
  labs(x="Numerator",
       y="Denominator",
       title = "Muscle Hypertrophy Model Comparision Matrix (Response Ratio)",
       subtitle = "2 x Log(Bayes Factor) Positive Values Favor Numerator",
       fill="2 x Log(BF)",
       caption = "Kass and Raferty (1995) Scale: -Inf to 0 = Negative; 0 to 2 = Weak; 2 to 6 = Positive; 6 to 10 = Strong; 10 to +Inf = Very Strong")+
  geom_tile(aes(fill=value))+
  geom_text(aes(label=value),color="black")+
  scale_fill_gradient2(low = "gray50",
                       high = "#0000a7",
                       mid="white",
                       midpoint = 0)+
  scale_x_discrete(labels=c("Linear\n(Intercept)",
                            "Linear\n(Slope)",
                            "Linear\n(Slopes)",
                            "Linear-log\n(Intercept)",
                            "Linear-log\n(Slope)",
                            "Linear-log\n(Slopes)",
                            "Quadratic\n(Intercept)",
                            "Quadratic\n(Slope)",
                            "Quadratic\n(Slopes)",
                            "Cubic\n(Intercept)",
                            "Cubic\n(Slope)",
                            "Cubic\n(Slopes)",
                            "Linear\nSpline\n(Intercept)",
                            "Linear\nSpline\n(Slope)",
                            "Linear\nSpline\n(Slopes)",
                            "Restricted\nCubic Spline\n(Intercept)",
                            "Restricted\nCubic Spline\n(Slope)",
                            "Restricted\nCubic Spline\n(Slopes)"))+
  scale_y_discrete(labels=rev(c("Linear\n(Intercept)",
                                "Linear\n(Slope)",
                                "Linear\n(Slopes)",
                                "Linear-log\n(Intercept)",
                                "Linear-log\n(Slope)",
                                "Linear-log\n(Slopes)",
                                "Quadratic\n(Intercept)",
                                "Quadratic\n(Slope)",
                                "Quadratic\n(Slopes)",
                                "Cubic\n(Intercept)",
                                "Cubic\n(Slope)",
                                "Cubic\n(Slopes)",
                                "Linear\nSpline\n(Intercept)",
                                "Linear\nSpline\n(Slope)",
                                "Linear\nSpline\n(Slopes)",
                                "Restricted\nCubic Spline\n(Intercept)",
                                "Restricted\nCubic Spline\n(Slope)", 
                                "Restricted\nCubic Spline\n(Slopes)"))) # linear.ri best fit

ggsave(primary.model.compare.hyp.rom,device = "pdf",filename="primary.model.compare.hyp.rom.pdf",width = 14,height = 9,
       path = "/Users/zacrobinson/Documents/For Publication PTF Meta/Plots")


#### refit rom model with robust variance estimates
primary.hyp.rom.best.model.prep=rma.mv(yi.rom,vi.rom,data = primary.data.hyp,
                                mods = ~ avg.rir+load.set+set.rep.equated+weeks+train.status,
                                random = list(~1|study,~1|group,~1|obs),
                                method = "REML", test = "t",dfs = "contain")

primary.hyp.rom.best.model=robust.rma.mv(primary.hyp.rom.best.model.prep,
                                       cluster = primary.data.hyp$study)

#### I2 and R2 for rom model
primary.hyp.rom.best.model.r2=r2_ml(primary.hyp.rom.best.model)*100
save(primary.hyp.rom.best.model.r2,file = "primary.hyp.rom.best.model.r2.RData")
primary.hyp.rom.best.model.i2=i2_ml(primary.hyp.rom.best.model)
save(primary.hyp.rom.best.model.i2,file = "primary.hyp.rom.best.model.i2.RData")

#### rom model output
primary.hyp.rom.best.output<-subset(as.data.frame(coef(summary(primary.hyp.rom.best.model))),select = c(-pval,-tval))
primary.hyp.rom.best.output.prep1<-qt(.975,primary.hyp.rom.best.output$df)
primary.hyp.rom.best.output.prep2<-sqrt((primary.hyp.rom.best.output$se^2)+sum(primary.hyp.rom.best.model$sigma2))
primary.hyp.rom.best.output$lower.PI<-primary.hyp.rom.best.output$estimate-(primary.hyp.rom.best.output.prep1*primary.hyp.rom.best.output.prep2)
primary.hyp.rom.best.output$upper.PI<-primary.hyp.rom.best.output$estimate+(primary.hyp.rom.best.output.prep1*primary.hyp.rom.best.output.prep2)
primary.hyp.rom.best.output$Parameter<-c("Intercept",
                                         "Estimated RIR",
                                         "Load",
                                         "Volume Equated [Rep]",
                                         "Volume Equated [Both]",
                                         "Weeks",
                                         "Training Status [Untrained]")

primary.hyp.rom.best.output.table<-relocate(primary.hyp.rom.best.output,Parameter)%>%
  gt()%>%
  cols_label(ci.lb="Lower Bound",
             ci.ub="Upper Bound",
             lower.PI="Lower Bound",
             upper.PI="Upper Bound",
             se="Standard Error",
             estimate="Estimate")%>%
  fmt_markdown(columns = 1)%>%
  tab_header(title = "Coefficients of Linear Multilevel Model for Muscle Hypertrophy (Response Ratio)")%>%
  tab_spanner(columns = 5:6,label = "95% Confidence Interval")%>%
  tab_spanner(columns = 7:8,label = "95% Predicition Interval")%>%
  opt_stylize(style = 1,color = "gray")%>%
  tab_options(column_labels.background.color ="gray97",column_labels.font.weight = "bold")%>%
  tab_style(style = cell_text(style = "italic"),locations = cells_row_groups())%>%
  tab_style(style = cell_text(align = "center"),locations = cells_column_labels())%>%
  tab_style(style = cell_text(weight  = "bold"),locations = cells_title())%>%
  tab_style(style = cell_text(align = "center"),locations = cells_body(columns=-1))%>%
  cols_width(-df~px(133))%>%
  cols_width(Parameter~px(200))

gtsave(primary.hyp.rom.best.output.table,filename = "primary.hyp.rom.best.output.table.png",
       path = "/Users/zacrobinson/Documents/For Publication PTF Meta/Models")

# MAIN EFFECTS ------------------------------------------------------

### strength
#### g marginal means
primary.str.g.best.emmeans=emmprep(primary.str.g.best.model,
                                     at=list(avg.rir=0:23))%>%
  emmeans(~avg.rir,
          weights="prop")%>%
  pred_interval_esmeans(model=primary.str.g.best.model)%>%
  as.data.frame()

#### g marginal slope
primary.str.g.best.slope.prep<-emmprep(primary.str.g.best.model,
                                       at=list(avg.rir=primary.str.mean.rir))

primary.str.g.best.slope<-confint(emmeans(primary.str.g.best.slope.prep,consec~avg.rir,
                                          weights = "prop")$contrasts)

primary.str.g.best.slope.prep2<-qt(.975,primary.str.g.best.slope$df)
primary.str.g.best.slope.prep3<-sqrt((primary.str.g.best.slope$SE^2)+sum(primary.str.g.best.model$sigma2))

primary.str.g.best.slope$lower.PI<-primary.str.g.best.slope$estimate-(primary.str.g.best.slope.prep2*primary.str.g.best.slope.prep3)
primary.str.g.best.slope$upper.PI<-primary.str.g.best.slope$estimate+(primary.str.g.best.slope.prep2*primary.str.g.best.slope.prep3)

save(primary.str.g.best.slope,file = "primary.str.g.best.slope.RData")

#### g plot
primary.str.g.best.model.plot=ggplot(data=primary.str.g.best.emmeans,
       aes(x=avg.rir,
           y=emmean))+THEME+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(data = primary.data.str,aes(x=avg.rir,y=yi,size=weights),alpha=0.75,shape=16,color="#c1272d")+
  geom_ribbon(aes(x=avg.rir,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(aes(x=avg.rir,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(linewidth=1,color="#c1272d")+
  scale_x_continuous(breaks = c(0,2,4,6,8,10,12,14,16,18,20,22))+
  scale_y_continuous(n.breaks = 9)+
  scale_size_continuous(range = c(2,8))+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Marginal Effects from Linear-log Multilevel Model for Maximal Strength (Standardized Mean Change)")

ggsave(primary.str.g.best.model.plot,device = "pdf",filename="primary.str.g.best.model.plot.pdf",width = 12,height = 8,
       path = "/Users/zacrobinson/Documents/For Publication PTF Meta/Plots")

#### g table
primary.str.g.best.model.table=primary.str.g.best.emmeans%>%
  gt()%>%
  tab_header(title = "Marginal Means from Linear-log Multilevel Model for Maximal Strength (Standardized Mean Change)")%>%
  cols_label(avg.rir="Estimated Proximity to Failure (RIR)",
             emmean="Effect Size (Hedge's g)",
             SE="Standard Error",
             lower.CL="Lower Bound",
             upper.CL="Upper Bound",
             lower.PI="Lower Bound",
             upper.PI="Upper Bound")%>%
  tab_spanner(label = "95% Confidence Interval",columns = 5:6)%>%
  tab_spanner(label = "95% Predicition Interval",columns = 7:8)%>%
  opt_stylize(style = 1,color = "gray")%>%
  tab_options(column_labels.background.color ="gray97",column_labels.font.weight = "bold")%>%
  tab_style(style = cell_text(style = "italic"),locations = cells_row_groups())%>%
  tab_style(style = cell_text(align = "center"),locations = cells_column_labels())%>%
  tab_style(style = cell_text(weight  = "bold"),locations = cells_title())%>%
  cols_width(-df~px(150))%>%
  tab_style(style = cell_text(align = "center"),locations = cells_body())

gtsave(primary.str.g.best.model.table,filename = "primary.str.g.best.model.table.png",vheight=1000,vwidth=1200,
       path = "/Users/zacrobinson/Documents/For Publication PTF Meta/Tables")

#### rom marginal means
primary.str.rom.best.emmeans=emmprep(primary.str.rom.best.model,
                                   at=list(avg.rir=0:23))%>%
  emmeans(~avg.rir,
          weights="prop")%>%
  pred_interval_esmeans(model=primary.str.rom.best.model)%>%
  as.data.frame()%>%
  mutate_at(c(2,3,5,6,7,8),expo)

#### rom marginal slope
primary.str.rom.best.slope.prep<-emmprep(primary.str.rom.best.model,
                                       at=list(avg.rir=primary.str.mean.rir))

primary.str.rom.best.slope<-confint(emmeans(primary.str.rom.best.slope.prep,consec~avg.rir,
                                          weights = "prop")$contrasts)

primary.str.rom.best.slope.prep2<-qt(.975,primary.str.rom.best.slope$df)
primary.str.rom.best.slope.prep3<-sqrt((primary.str.rom.best.slope$SE^2)+sum(primary.str.rom.best.model$sigma2))

primary.str.rom.best.slope$lower.PI<-primary.str.rom.best.slope$estimate-(primary.str.rom.best.slope.prep2*primary.str.rom.best.slope.prep3)
primary.str.rom.best.slope$upper.PI<-primary.str.rom.best.slope$estimate+(primary.str.rom.best.slope.prep2*primary.str.rom.best.slope.prep3)

primary.str.rom.best.slope=primary.str.rom.best.slope%>%
  mutate_at(c(2,3,5,6,7,8),expo)

save(primary.str.rom.best.slope,file = "primary.str.rom.best.slope.RData")

#### rom plot
primary.str.rom.best.model.plot=ggplot(data=primary.str.rom.best.emmeans,
       aes(x=avg.rir,
           y=emmean))+THEME+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(data = primary.data.str,aes(x=avg.rir,y=yi.rom.exp,size=weights.rom),alpha=0.75,shape=16,color="#c1272d")+
  geom_ribbon(aes(x=avg.rir,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(aes(x=avg.rir,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(linewidth=1,color="#c1272d")+
  scale_x_continuous(breaks = c(0,2,4,6,8,10,12,14,16,18,20,22))+
  scale_y_continuous(n.breaks = 9)+
  scale_size_continuous(range = c(2,8))+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio for Muscle Size (% Change)",
       title = "Marginal Effect from Linear Multilevel Model for Maximal Strength (Response Ratio)")

ggsave(primary.str.rom.best.model.plot,device = "pdf",filename="primary.str.rom.best.model.plot.pdf",width = 12,height = 8,
       path = "/Users/zacrobinson/Documents/For Publication PTF Meta/Plots")

#### rom table
primary.str.rom.best.model.table=primary.str.rom.best.emmeans%>%
  gt()%>%
  tab_header(title = "Marginal Means from Linear Multilevel Model for Maximal Strength (Response Ratio)")%>%
  cols_label(avg.rir="Estimated Proximity to Failure (RIR)",
             emmean="Effect Size (% Change)",
             SE="Standard Error",
             lower.CL="Lower Bound",
             upper.CL="Upper Bound",
             lower.PI="Lower Bound",
             upper.PI="Upper Bound")%>%
  tab_spanner(label = "95% Confidence Interval",columns = 5:6)%>%
  tab_spanner(label = "95% Predicition Interval",columns = 7:8)%>%
  opt_stylize(style = 1,color = "gray")%>%
  tab_options(column_labels.background.color ="gray97",column_labels.font.weight = "bold")%>%
  tab_style(style = cell_text(style = "italic"),locations = cells_row_groups())%>%
  tab_style(style = cell_text(align = "center"),locations = cells_column_labels())%>%
  tab_style(style = cell_text(weight  = "bold"),locations = cells_title())%>%
  cols_width(-df~px(150))%>%
  tab_style(style = cell_text(align = "center"),locations = cells_body())

gtsave(primary.str.rom.best.model.table,filename = "primary.str.rom.best.model.table.png",vheight=1000,vwidth=1200,
       path = "/Users/zacrobinson/Documents/For Publication PTF Meta/Tables")

### hypertrophy
#### g marginal means
primary.hyp.g.best.emmeans=emmprep(primary.hyp.g.best.model,
                                   at=list(avg.rir=0:23))%>%
  emmeans(~avg.rir,
          weights="prop")%>%
  pred_interval_esmeans(model=primary.hyp.g.best.model)%>%
  as.data.frame()

#### g marginal slope
primary.hyp.g.best.slope.prep<-emmprep(primary.hyp.g.best.model,
                                       at=list(avg.rir=primary.hyp.mean.rir))

primary.hyp.g.best.slope<-confint(emmeans(primary.hyp.g.best.slope.prep,consec~avg.rir,
                                          weights = "prop")$contrasts)

primary.hyp.g.best.slope.prep2<-qt(.975,primary.hyp.g.best.slope$df)
primary.hyp.g.best.slope.prep3<-sqrt((primary.hyp.g.best.slope$SE^2)+sum(primary.hyp.g.best.model$sigma2))

primary.hyp.g.best.slope$lower.PI<-primary.hyp.g.best.slope$estimate-(primary.hyp.g.best.slope.prep2*primary.hyp.g.best.slope.prep3)
primary.hyp.g.best.slope$upper.PI<-primary.hyp.g.best.slope$estimate+(primary.hyp.g.best.slope.prep2*primary.hyp.g.best.slope.prep3)

save(primary.hyp.g.best.slope,file = "primary.hyp.g.best.slope.RData")

#### g plot
primary.hyp.g.best.model.plot=ggplot(data=primary.hyp.g.best.emmeans,
       aes(x=avg.rir,
           y=emmean))+THEME+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(data = primary.data.hyp,aes(x=avg.rir,y=yi,size=weights),alpha=0.75,shape=16,color="#0000a7")+
  geom_ribbon(aes(x=avg.rir,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(aes(x=avg.rir,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(linewidth=1,color="#0000a7")+
  scale_x_continuous(breaks = c(0,2,4,6,8,10,12,14,16,18,20,22))+
  scale_y_continuous(n.breaks = 9)+
  scale_size_continuous(range = c(2,8))+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Muscle Size (Hedge's g)",
       title = "Marginal Effects from Linear Multilevel Model for Muscle Hypertrophy (Standardized Mean Change)")

ggsave(primary.hyp.g.best.model.plot,device = "pdf",filename="primary.hyp.g.best.model.plot.pdf",width = 12,height = 8,
       path = "/Users/zacrobinson/Documents/For Publication PTF Meta/Plots")

#### g table
primary.hyp.g.best.model.table=primary.hyp.g.best.emmeans%>%
  gt()%>%
  tab_header(title = "Marginal Means from Linear Multilevel Model for Muscle Hypertrophy (Standardized Mean Change)")%>%
  cols_label(avg.rir="Estimated Proximity to Failure (RIR)",
             emmean="Effect Size (Hedge's g)",
             SE="Standard Error",
             lower.CL="Lower Bound",
             upper.CL="Upper Bound",
             lower.PI="Lower Bound",
             upper.PI="Upper Bound")%>%
  tab_spanner(label = "95% Confidence Interval",columns = 5:6)%>%
  tab_spanner(label = "95% Predicition Interval",columns = 7:8)%>%
  opt_stylize(style = 1,color = "gray")%>%
  tab_options(column_labels.background.color ="gray97",column_labels.font.weight = "bold")%>%
  tab_style(style = cell_text(style = "italic"),locations = cells_row_groups())%>%
  tab_style(style = cell_text(align = "center"),locations = cells_column_labels())%>%
  tab_style(style = cell_text(weight  = "bold"),locations = cells_title())%>%
  cols_width(-df~px(150))%>%
  tab_style(style = cell_text(align = "center"),locations = cells_body())

gtsave(primary.hyp.g.best.model.table,filename = "primary.hyp.g.best.model.table.png",vheight=1000,vwidth=1200,
       path = "/Users/zacrobinson/Documents/For Publication PTF Meta/Tables")



#### rom marginal means
primary.hyp.rom.best.emmeans=emmprep(primary.hyp.rom.best.model,
                                     at=list(avg.rir=0:23))%>%
  emmeans(~avg.rir,
          weights="prop")%>%
  pred_interval_esmeans(model=primary.hyp.rom.best.model)%>%
  as.data.frame()%>%
  mutate_at(c(2,3,5,6,7,8),expo)

#### rom marginal slope
primary.hyp.rom.best.slope.prep<-emmprep(primary.hyp.rom.best.model,
                                         at=list(avg.rir=primary.hyp.mean.rir))

primary.hyp.rom.best.slope<-confint(emmeans(primary.hyp.rom.best.slope.prep,consec~avg.rir,
                                            weights = "prop")$contrasts)

primary.hyp.rom.best.slope.prep2<-qt(.975,primary.hyp.rom.best.slope$df)
primary.hyp.rom.best.slope.prep3<-sqrt((primary.hyp.rom.best.slope$SE^2)+sum(primary.hyp.rom.best.model$sigma2))

primary.hyp.rom.best.slope$lower.PI<-primary.hyp.rom.best.slope$estimate-(primary.hyp.rom.best.slope.prep2*primary.hyp.rom.best.slope.prep3)
primary.hyp.rom.best.slope$upper.PI<-primary.hyp.rom.best.slope$estimate+(primary.hyp.rom.best.slope.prep2*primary.hyp.rom.best.slope.prep3)

primary.hyp.rom.best.slope=primary.hyp.rom.best.slope%>%
  mutate_at(c(2,3,5,6,7,8),expo)

save(primary.hyp.rom.best.slope,file = "primary.hyp.rom.best.slope.RData")

#### rom plot
primary.hyp.rom.best.model.plot=ggplot(data=primary.hyp.rom.best.emmeans,
       aes(x=avg.rir,
           y=emmean))+THEME+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(data = primary.data.hyp,aes(x=avg.rir,y=yi.rom.exp,size=weights.rom),alpha=0.75,shape=16,color="#0000a7")+
  geom_ribbon(aes(x=avg.rir,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(aes(x=avg.rir,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(linewidth=1,color="#0000a7")+
  scale_x_continuous(breaks = c(0,2,4,6,8,10,12,14,16,18,20,22))+
  scale_y_continuous(n.breaks = 9)+
  scale_size_continuous(range = c(2,8))+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio for Muscle Size (% Change)",
       title = "Marginal Effects from Linear Multilevel Model for Muscle Hypertrophy (Response Ratio)")

ggsave(primary.hyp.rom.best.model.plot,device = "pdf",filename="primary.hyp.rom.best.model.plot.pdf",width = 12,height = 8,
       path = "/Users/zacrobinson/Documents/For Publication PTF Meta/Plots")

#### rom table
primary.hyp.rom.best.model.table=primary.hyp.rom.best.emmeans%>%
  gt()%>%
  tab_header(title = "Marginal Means from Linear Multilevel Model for Muscle Hypertrophy (Response Ratio)")%>%
  cols_label(avg.rir="Estimated Proximity to Failure (RIR)",
             emmean="Effect Size (% Change)",
             SE="Standard Error",
             lower.CL="Lower Bound",
             upper.CL="Upper Bound",
             lower.PI="Lower Bound",
             upper.PI="Upper Bound")%>%
  tab_spanner(label = "95% Confidence Interval",columns = 5:6)%>%
  tab_spanner(label = "95% Predicition Interval",columns = 7:8)%>%
  opt_stylize(style = 1,color = "gray")%>%
  tab_options(column_labels.background.color ="gray97",column_labels.font.weight = "bold")%>%
  tab_style(style = cell_text(style = "italic"),locations = cells_row_groups())%>%
  tab_style(style = cell_text(align = "center"),locations = cells_column_labels())%>%
  tab_style(style = cell_text(weight  = "bold"),locations = cells_title())%>%
  cols_width(-df~px(150))%>%
  tab_style(style = cell_text(align = "center"),locations = cells_body())

gtsave(primary.hyp.rom.best.model.table,filename = "primary.hyp.rom.best.model.table.png",vheight=1000,vwidth=1200,
       path = "/Users/zacrobinson/Documents/For Publication PTF Meta/Tables")

# MODERATORS --------------------------------------------------------------

### strength
#### g load
primary.str.g.moderator.load=primary.str.g.best.model.prep%>%
  update(~. + log1p(avg.rir):load.set)%>%
  robust.rma.mv(cluster = primary.data.str$study)
  
primary.str.g.moderator.load.emms=primary.str.g.moderator.load%>%
  emmprep(at=list(avg.rir=0:23,
                  load.set=c(30,60,90)))%>%
  emmeans(~avg.rir+load.set,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.str.g.moderator.load)

primary.str.g.moderator.load.contrast=emmprep(primary.str.g.moderator.load,
                                              at=list(avg.rir=primary.str.mean.rir,
                                                      load.set=primary.str.mean.load.set))%>%
  emmeans(~avg.rir+load.set,
          weights = "prop")%>%
  contrast("consec",interaction=TRUE)%>%
  confint()

primary.str.g.moderator.load.contrast.label=paste0('Interaction Contrast (Load) = ',
                                                   bquote(.(round(primary.str.g.moderator.load.contrast[[3]],5))),
                                                   " [95% CI: ",bquote(.(round(primary.str.g.moderator.load.contrast[[6]],5))),
                                                   ", ",bquote(.(round(primary.str.g.moderator.load.contrast[[7]],5))),"]")
  
primary.str.g.moderator.load.plot=ggplot(data = primary.data.str%>%
         mutate(load.set=ifelse(load.set<45,"30% of 1RM",
                                ifelse(load.set>45&load.set<75,"60% of 1RM","90% of 1RM"))),
       aes(x=avg.rir,
           y=yi))+THEME+theme(plot.caption = element_text(face = "italic"))+
  facet_wrap(~load.set)+
  labs(x="Estimated Proximity to Failure (RIR)",
                            y="Standardized Mean Change in Strength (Hedge's g)",
                            title = "Load (% of 1RM) x Estimated RIR",
       caption = primary.str.g.moderator.load.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = primary.str.g.moderator.load.emms%>%
                mutate(load.set=ifelse(load.set==30,"30% of 1RM",
                                       ifelse(load.set==60,"60% of 1RM","90% of 1RM"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = primary.str.g.moderator.load.emms%>%
                mutate(load.set=ifelse(load.set==30,"30% of 1RM",
                                       ifelse(load.set==60,"60% of 1RM","90% of 1RM"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = primary.str.g.moderator.load.emms%>%
              mutate(load.set=ifelse(load.set==30,"30% of 1RM",
                                     ifelse(load.set==60,"60% of 1RM","90% of 1RM"))),
    aes(x=avg.rir,
        y=emmean),
    linewidth=1,
    color="#c1272d")+
  scale_size_continuous(guide="none")

#### g volume equating method
primary.str.g.moderator.volume.method=primary.str.g.best.model.prep%>%
  update(~. + log1p(avg.rir):set.rep.equated)%>%
  robust.rma.mv(cluster = primary.data.str$study)

primary.str.g.moderator.volume.method.emms=primary.str.g.moderator.volume.method%>%emmprep(at=list(avg.rir=0:23))%>%
  emmeans(~avg.rir+set.rep.equated,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.str.g.moderator.volume.method)

primary.str.g.moderator.volume.method.contrast=emmprep(primary.str.g.moderator.volume.method,
                                              at=list(avg.rir=primary.str.mean.rir))%>%
  emmeans(~avg.rir+set.rep.equated,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()

primary.str.g.moderator.volume.method.contrast.label=paste0('Interaction Contrast (Set - Rep) = ',
                                                   bquote(.(round(primary.str.g.moderator.volume.method.contrast[[1,3]],5))),
                                                   " [95% CI: ",bquote(.(round(primary.str.g.moderator.volume.method.contrast[[1,6]],5))),
                                                   ", ",bquote(.(round(primary.str.g.moderator.volume.method.contrast[[1,7]],5))),"]\n",
                                                   'Interaction Contrast (Set - Both) = ',
                                                   bquote(.(round(primary.str.g.moderator.volume.method.contrast[[2,3]],5))),
                                                   " [95% CI: ",bquote(.(round(primary.str.g.moderator.volume.method.contrast[[2,6]],5))),
                                                   ", ",bquote(.(round(primary.str.g.moderator.volume.method.contrast[[2,7]],5))),"]\n",
                                                   'Interaction Contrast (Rep - Both) = ',
                                                   bquote(.(round(primary.str.g.moderator.volume.method.contrast[[3,3]],5))),
                                                   " [95% CI: ",bquote(.(round(primary.str.g.moderator.volume.method.contrast[[3,6]],5))),
                                                   ", ",bquote(.(round(primary.str.g.moderator.volume.method.contrast[[3,7]],5))),"]")


primary.str.g.moderator.volume.method.plot=ggplot(data = primary.data.str%>%
         mutate(set.rep.equated=ifelse(set.rep.equated=="set","Set Equated",
                                       ifelse(set.rep.equated=="rep","Repetition Equated","Both Equated"))),
       aes(x=avg.rir,
           y=yi))+THEME+theme(legend.position = "right",
                              plot.caption = element_text(face = "italic"))+
  facet_wrap(~set.rep.equated)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Method of Volume Equating x Estimated RIR",
       caption = primary.str.g.moderator.volume.method.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = primary.str.g.moderator.volume.method.emms%>%
                mutate(set.rep.equated=ifelse(set.rep.equated=="set","Set Equated",
                                              ifelse(set.rep.equated=="rep","Repetition Equated","Both Equated"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = primary.str.g.moderator.volume.method.emms%>%
                mutate(set.rep.equated=ifelse(set.rep.equated=="set","Set Equated",
                                              ifelse(set.rep.equated=="rep","Repetition Equated","Both Equated"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = primary.str.g.moderator.volume.method.emms%>%
              mutate(set.rep.equated=ifelse(set.rep.equated=="set","Set Equated",
                                            ifelse(set.rep.equated=="rep","Repetition Equated","Both Equated"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### g training status
primary.str.g.moderator.train.status=primary.str.g.best.model.prep%>%
  update(~. + log1p(avg.rir):train.status)%>%
  robust.rma.mv(cluster = primary.data.str$study)

primary.str.g.moderator.train.status.emms=primary.str.g.moderator.train.status%>%emmprep(at=list(avg.rir=0:23))%>%
  emmeans(~avg.rir+train.status,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.str.g.moderator.train.status)

primary.str.g.moderator.train.status.contrast=emmprep(primary.str.g.moderator.train.status,
                                                       at=list(avg.rir=primary.str.mean.rir))%>%
  emmeans(~avg.rir+train.status,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()

primary.str.g.moderator.train.status.contrast.label=paste0('Interaction Contrast (Trained - Untrained) = ',
                                                            bquote(.(round(primary.str.g.moderator.train.status.contrast[[1,3]],5))),
                                                            " [95% CI: ",bquote(.(round(primary.str.g.moderator.train.status.contrast[[1,6]],5))),
                                                            ", ",bquote(.(round(primary.str.g.moderator.train.status.contrast[[1,7]],5))),"]")


primary.str.g.moderator.train.status.plot=ggplot(data = primary.data.str%>%
                                                    mutate(train.status=ifelse(train.status=="trained","Trained","Untrained")),
                                                  aes(x=avg.rir,
                                                      y=yi))+THEME+theme(legend.position = "right",
                                                                         plot.caption = element_text(face = "italic"))+
  facet_wrap(~train.status)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Training Status x Estimated RIR",
       caption = primary.str.g.moderator.train.status.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = primary.str.g.moderator.train.status.emms%>%
                mutate(train.status=ifelse(train.status=="trained","Trained","Untrained")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = primary.str.g.moderator.train.status.emms%>%
                mutate(train.status=ifelse(train.status=="trained","Trained","Untrained")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = primary.str.g.moderator.train.status.emms%>%
              mutate(train.status=ifelse(train.status=="trained","Trained","Untrained")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### g weeks
primary.str.g.moderator.weeks=primary.str.g.best.model.prep%>%
  update(~. + log1p(avg.rir):weeks)%>%
  robust.rma.mv(cluster = primary.data.str$study)

primary.str.g.moderator.weeks.emms=primary.str.g.moderator.weeks%>%
  emmprep(at=list(avg.rir=0:23,
                  weeks=c(4,8,12)))%>%
  emmeans(~avg.rir+weeks,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.str.g.moderator.weeks)

primary.str.g.moderator.weeks.contrast=emmprep(primary.str.g.moderator.weeks,
                                              at=list(avg.rir=primary.str.mean.rir,
                                                      weeks=primary.str.mean.weeks))%>%
  emmeans(~avg.rir+weeks,
          weights = "prop")%>%
  contrast("consec",interaction=TRUE)%>%
  confint()

primary.str.g.moderator.weeks.contrast.label=paste0('Interaction Contrast (Weeks) = ',
                                                   bquote(.(round(primary.str.g.moderator.weeks.contrast[[3]],5))),
                                                   " [95% CI: ",bquote(.(round(primary.str.g.moderator.weeks.contrast[[6]],5))),
                                                   ", ",bquote(.(round(primary.str.g.moderator.weeks.contrast[[7]],5))),"]")

primary.str.g.moderator.weeks.plot=ggplot(data = primary.data.str%>%
                                           mutate(weeks=ifelse(weeks<6,"4 Weeks",
                                                                  ifelse(weeks>6&weeks<10,"8 Weeks","12 Weeks"))),
                                         aes(x=avg.rir,
                                             y=yi))+THEME+theme(plot.caption = element_text(face = "italic"))+
  facet_wrap(~fct_relevel(weeks,c("4 Weeks","8 Weeks","12 Weeks")))+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Intervention Duration (Weeks) x Estimated RIR",
       caption = primary.str.g.moderator.weeks.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = primary.str.g.moderator.weeks.emms%>%
                mutate(weeks=ifelse(weeks==4,"4 Weeks",
                                       ifelse(weeks==8,"8 Weeks","12 Weeks"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = primary.str.g.moderator.weeks.emms%>%
                mutate(weeks=ifelse(weeks==4,"4 Weeks",
                                    ifelse(weeks==8,"8 Weeks","12 Weeks"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = primary.str.g.moderator.weeks.emms%>%
              mutate(weeks=ifelse(weeks==4,"4 Weeks",
                                  ifelse(weeks==8,"8 Weeks","12 Weeks"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### g adjusted sets 
primary.str.g.moderator.adj.sets=primary.str.g.best.model.prep%>%
  update(~. + log1p(avg.rir)*adj.sets.week)%>%
  robust.rma.mv(cluster = primary.data.str$study)

primary.str.g.moderator.adj.sets.emms=primary.str.g.moderator.adj.sets%>%
  emmprep(at=list(avg.rir=0:23,
                  adj.sets.week=c(8,16,24)))%>%
  emmeans(~avg.rir+adj.sets.week,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.str.g.moderator.adj.sets)

primary.str.g.moderator.adj.sets.contrast=emmprep(primary.str.g.moderator.adj.sets,
                                               at=list(avg.rir=primary.str.mean.rir,
                                                       adj.sets.week=primary.str.mean.adj.sets.week))%>%
  emmeans(~avg.rir+adj.sets.week,
          weights = "prop")%>%
  contrast("consec",interaction=TRUE)%>%
  confint()

primary.str.g.moderator.adj.sets.contrast.label=paste0('Interaction Contrast (Adjusted Direct Sets) = ',
                                                    bquote(.(round(primary.str.g.moderator.adj.sets.contrast[[3]],5))),
                                                    " [95% CI: ",bquote(.(round(primary.str.g.moderator.adj.sets.contrast[[6]],5))),
                                                    ", ",bquote(.(round(primary.str.g.moderator.adj.sets.contrast[[7]],5))),"]")

primary.str.g.moderator.adj.sets.plot=ggplot(data = primary.data.str%>%
                                            mutate(adj.sets.week=ifelse(adj.sets.week<12,"8 Sets",
                                                                ifelse(adj.sets.week>12&adj.sets.week<20,"16 Sets","24 Sets"))),
                                          aes(x=avg.rir,
                                              y=yi))+THEME+theme(plot.caption = element_text(face = "italic"))+
  facet_wrap(~fct_relevel(adj.sets.week,c("8 Sets","16 Sets","24 Sets")))+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Adjusted Direct Sets Per Week Per Muscle Group x Estimated RIR",
       caption = primary.str.g.moderator.adj.sets.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = primary.str.g.moderator.adj.sets.emms%>%
                mutate(adj.sets.week=ifelse(adj.sets.week==8,"8 Sets",
                                       ifelse(adj.sets.week==16,"16 Sets","24 Sets"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = primary.str.g.moderator.adj.sets.emms%>%
                mutate(adj.sets.week=ifelse(adj.sets.week==8,"8 Sets",
                                       ifelse(adj.sets.week==16,"16 Sets","24 Sets"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = primary.str.g.moderator.adj.sets.emms%>%
              mutate(adj.sets.week=ifelse(adj.sets.week==8,"8 Sets",
                                     ifelse(adj.sets.week==16,"16 Sets","24 Sets"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### g frequency
primary.str.g.moderator.frequency.per.muscle=primary.str.g.best.model.prep%>%
  update(~. + log1p(avg.rir)*frequency.per.muscle)%>%
  robust.rma.mv(cluster = primary.data.str$study)

primary.str.g.moderator.frequency.per.muscle.emms=primary.str.g.moderator.frequency.per.muscle%>%
  emmprep(at=list(avg.rir=0:23,
                  frequency.per.muscle=c(1,2,3)))%>%
  emmeans(~avg.rir+frequency.per.muscle,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.str.g.moderator.frequency.per.muscle)

primary.str.g.moderator.frequency.per.muscle.contrast=emmprep(primary.str.g.moderator.frequency.per.muscle,
                                                  at=list(avg.rir=primary.str.mean.rir,
                                                          frequency.per.muscle=primary.str.mean.frequency.per.muscle))%>%
  emmeans(~avg.rir+frequency.per.muscle,
          weights = "prop")%>%
  contrast("consec",interaction=TRUE)%>%
  confint()

primary.str.g.moderator.frequency.per.muscle.contrast.label=paste0('Interaction Contrast (Direct Frequency) = ',
                                                       bquote(.(round(primary.str.g.moderator.frequency.per.muscle.contrast[[3]],5))),
                                                       " [95% CI: ",bquote(.(round(primary.str.g.moderator.frequency.per.muscle.contrast[[6]],5))),
                                                       ", ",bquote(.(round(primary.str.g.moderator.frequency.per.muscle.contrast[[7]],5))),"]")

primary.str.g.moderator.frequency.per.muscle.plot=ggplot(data = primary.data.str%>%
                                               mutate(frequency.per.muscle=ifelse(frequency.per.muscle<1.5,"1x Per Week",
                                                                           ifelse(frequency.per.muscle>1.5&frequency.per.muscle<2.5,"2x Per Week","3x Per Week"))),
                                             aes(x=avg.rir,
                                                 y=yi))+THEME+theme(plot.caption = element_text(face = "italic"))+
  facet_wrap(~frequency.per.muscle)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Direct Frequency Per Week Per Muscle Group x Estimated RIR",
       caption = primary.str.g.moderator.frequency.per.muscle.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = primary.str.g.moderator.frequency.per.muscle.emms%>%
                mutate(frequency.per.muscle=ifelse(frequency.per.muscle==1,"1x Per Week",
                                       ifelse(frequency.per.muscle==2,"2x Per Week","3x Per Week"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = primary.str.g.moderator.frequency.per.muscle.emms%>%
                mutate(frequency.per.muscle=ifelse(frequency.per.muscle==1,"1x Per Week",
                                       ifelse(frequency.per.muscle==2,"2x Per Week","3x Per Week"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = primary.str.g.moderator.frequency.per.muscle.emms%>%
              mutate(frequency.per.muscle=ifelse(frequency.per.muscle==1,"1x Per Week",
                                     ifelse(frequency.per.muscle==2,"2x Per Week","3x Per Week"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### g age
primary.str.g.moderator.age=primary.str.g.best.model.prep%>%
  update(~. + log1p(avg.rir)*age)%>%
  robust.rma.mv(cluster = primary.data.str$study)

primary.str.g.moderator.age.emms=primary.str.g.moderator.age%>%
  emmprep(at=list(avg.rir=0:23,
                  age=c(20,40,60)))%>%
  emmeans(~avg.rir+age,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.str.g.moderator.age)

primary.str.g.moderator.age.contrast=emmprep(primary.str.g.moderator.age,
                                                              at=list(avg.rir=primary.str.mean.rir,
                                                                      age=primary.str.mean.age))%>%
  emmeans(~avg.rir+age,
          weights = "prop")%>%
  contrast("consec",interaction=TRUE)%>%
  confint()

primary.str.g.moderator.age.contrast.label=paste0('Interaction Contrast (Age) = ',
                                                                   bquote(.(round(primary.str.g.moderator.age.contrast[[3]],5))),
                                                                   " [95% CI: ",bquote(.(round(primary.str.g.moderator.age.contrast[[6]],5))),
                                                                   ", ",bquote(.(round(primary.str.g.moderator.age.contrast[[7]],5))),"]")

primary.str.g.moderator.age.plot=ggplot(data = primary.data.str%>%
                                                           mutate(age=ifelse(age<30,"20 Years Old",
                                                                                              ifelse(age>30&age<50,"40 Years Old","60 Years Old"))),
                                                         aes(x=avg.rir,
                                                             y=yi))+THEME+theme(plot.caption = element_text(face = "italic"))+
  facet_wrap(~age)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Age x Estimated RIR",
       caption = primary.str.g.moderator.age.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = primary.str.g.moderator.age.emms%>%
                mutate(age=ifelse(age==20,"20 Years Old",
                                       ifelse(age==40,"40 Years Old","60 Years Old"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = primary.str.g.moderator.age.emms%>%
                mutate(age=ifelse(age==20,"20 Years Old",
                                       ifelse(age==40,"40 Years Old","60 Years Old"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = primary.str.g.moderator.age.emms%>%
              mutate(age=ifelse(age==20,"20 Years Old",
                                     ifelse(age==40,"40 Years Old","60 Years Old"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### g sex
primary.str.g.moderator.sex.percent.male=primary.str.g.best.model.prep%>%
  update(~. + log1p(avg.rir)*sex.percent.male)%>%
  robust.rma.mv(cluster = primary.data.str$study)

primary.str.g.moderator.sex.percent.male.emms=primary.str.g.moderator.sex.percent.male%>%
  emmprep(at=list(avg.rir=0:23,
                  sex.percent.male=c(0,100)))%>%
  emmeans(~avg.rir+sex.percent.male,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.str.g.moderator.sex.percent.male)

primary.str.g.moderator.sex.percent.male.contrast=emmprep(primary.str.g.moderator.sex.percent.male,
                                             at=list(avg.rir=primary.str.mean.rir,
                                                     sex.percent.male=primary.str.mean.sex))%>%
  emmeans(~avg.rir+sex.percent.male,
          weights = "prop")%>%
  contrast("consec",interaction=TRUE)%>%
  confint()

primary.str.g.moderator.sex.percent.male.contrast.label=paste0('Interaction Contrast (Sex) = ',
                                                  bquote(.(round(primary.str.g.moderator.sex.percent.male.contrast[[3]],5))),
                                                  " [95% CI: ",bquote(.(round(primary.str.g.moderator.sex.percent.male.contrast[[6]],5))),
                                                  ", ",bquote(.(round(primary.str.g.moderator.sex.percent.male.contrast[[7]],5))),"]")

primary.str.g.moderator.sex.percent.male.plot=ggplot(data = primary.data.str%>%
                                          mutate(sex.percent.male=ifelse(sex.percent.male<50,"Female","Male")),
                                        aes(x=avg.rir,
                                            y=yi))+THEME+theme(plot.caption = element_text(face = "italic"))+
  facet_wrap(~sex.percent.male)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Sex (% Male) x Estimated RIR",
       caption = primary.str.g.moderator.sex.percent.male.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = primary.str.g.moderator.sex.percent.male.emms%>%
                mutate(sex.percent.male=ifelse(sex.percent.male==0,"Female","Male")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = primary.str.g.moderator.sex.percent.male.emms%>%
                mutate(sex.percent.male=ifelse(sex.percent.male==0,"Female","Male")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = primary.str.g.moderator.sex.percent.male.emms%>%
              mutate(sex.percent.male=ifelse(sex.percent.male==0,"Female","Male")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### g study design
primary.str.g.moderator.within.between.design=primary.str.g.best.model.prep%>%
  update(~. + log1p(avg.rir)*within.between.design)%>%
  robust.rma.mv(cluster = primary.data.str$study)

primary.str.g.moderator.within.between.design.emms=primary.str.g.moderator.within.between.design%>%emmprep(at=list(avg.rir=0:23))%>%
  emmeans(~avg.rir+within.between.design,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.str.g.moderator.within.between.design)

primary.str.g.moderator.within.between.design.contrast=emmprep(primary.str.g.moderator.within.between.design,
                                                      at=list(avg.rir=primary.str.mean.rir))%>%
  emmeans(~avg.rir+within.between.design,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()

primary.str.g.moderator.within.between.design.contrast.label=paste0('Interaction Contrast (Between - Within) = ',
                                                           bquote(.(round(primary.str.g.moderator.within.between.design.contrast[[1,3]],5))),
                                                           " [95% CI: ",bquote(.(round(primary.str.g.moderator.within.between.design.contrast[[1,6]],5))),
                                                           ", ",bquote(.(round(primary.str.g.moderator.within.between.design.contrast[[1,7]],5))),"]")


primary.str.g.moderator.within.between.design.plot=ggplot(data = primary.data.str%>%
                                                   mutate(within.between.design=ifelse(within.between.design=="Between","Between Subjects","Within Subjects")),
                                                 aes(x=avg.rir,
                                                     y=yi))+THEME+theme(legend.position = "right",
                                                                        plot.caption = element_text(face = "italic"))+
  facet_wrap(~within.between.design)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Study Design x Estimated RIR",
       caption = primary.str.g.moderator.within.between.design.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = primary.str.g.moderator.within.between.design.emms%>%
                mutate(within.between.design=ifelse(within.between.design=="Between","Between Subjects","Within Subjects")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = primary.str.g.moderator.within.between.design.emms%>%
                mutate(within.between.design=ifelse(within.between.design=="Between","Between Subjects","Within Subjects")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = primary.str.g.moderator.within.between.design.emms%>%
              mutate(within.between.design=ifelse(within.between.design=="Between","Between Subjects","Within Subjects")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### g upper / lower outcome
primary.str.g.moderator.upper.lower.other=primary.str.g.best.model.prep%>%
  update(~. + log1p(avg.rir)*upper.lower.other)%>%
  robust.rma.mv(cluster = primary.data.str$study)

primary.str.g.moderator.upper.lower.other.emms=primary.str.g.moderator.upper.lower.other%>%emmprep(at=list(avg.rir=0:23))%>%
  emmeans(~avg.rir+upper.lower.other,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.str.g.moderator.upper.lower.other)

primary.str.g.moderator.upper.lower.other.contrast=emmprep(primary.str.g.moderator.upper.lower.other,
                                                               at=list(avg.rir=primary.str.mean.rir))%>%
  emmeans(~avg.rir+upper.lower.other,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()

primary.str.g.moderator.upper.lower.other.contrast.label=paste0('Interaction Contrast (Lower - Upper) = ',
                                                                    bquote(.(round(primary.str.g.moderator.upper.lower.other.contrast[[1,3]],5))),
                                                                    " [95% CI: ",bquote(.(round(primary.str.g.moderator.upper.lower.other.contrast[[1,6]],5))),
                                                                    ", ",bquote(.(round(primary.str.g.moderator.upper.lower.other.contrast[[1,7]],5))),"]")


primary.str.g.moderator.upper.lower.other.plot=ggplot(data = primary.data.str%>%
                                                            mutate(upper.lower.other=ifelse(upper.lower.other=="Lower","Lower Body Outcomes","Upper Body Outcomes")),
                                                          aes(x=avg.rir,
                                                              y=yi))+THEME+theme(legend.position = "right",
                                                                                 plot.caption = element_text(face = "italic"))+
  facet_wrap(~upper.lower.other)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Upper/Lower Body Outcome x Estimated RIR",
       caption = primary.str.g.moderator.upper.lower.other.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = primary.str.g.moderator.upper.lower.other.emms%>%
                mutate(upper.lower.other=ifelse(upper.lower.other=="Lower","Lower Body Outcomes","Upper Body Outcomes")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = primary.str.g.moderator.upper.lower.other.emms%>%
                mutate(upper.lower.other=ifelse(upper.lower.other=="Lower","Lower Body Outcomes","Upper Body Outcomes")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = primary.str.g.moderator.upper.lower.other.emms%>%
              mutate(upper.lower.other=ifelse(upper.lower.other=="Lower","Lower Body Outcomes","Upper Body Outcomes")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### g exercise selection
primary.str.g.moderator.train.exercise=primary.str.g.best.model.prep%>%
  update(~. + log1p(avg.rir)*train.exercise)%>%
  robust.rma.mv(cluster = primary.data.str$study)

primary.str.g.moderator.train.exercise.emms=primary.str.g.moderator.train.exercise%>%emmprep(at=list(avg.rir=0:23))%>%
  emmeans(~avg.rir+train.exercise,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.str.g.moderator.train.exercise)

primary.str.g.moderator.train.exercise.contrast=emmprep(primary.str.g.moderator.train.exercise,
                                                           at=list(avg.rir=primary.str.mean.rir))%>%
  emmeans(~avg.rir+train.exercise,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()

primary.str.g.moderator.train.exercise.contrast.label=paste0('Interaction Contrast (Both - Multi-joint) = ',
                                                                bquote(.(round(primary.str.g.moderator.train.exercise.contrast[[1,3]],5))),
                                                                " [95% CI: ",bquote(.(round(primary.str.g.moderator.train.exercise.contrast[[1,6]],5))),
                                                                ", ",bquote(.(round(primary.str.g.moderator.train.exercise.contrast[[1,7]],5))),"]\n",
                                                                'Interaction Contrast (Both - Single-joint) = ',
                                                                bquote(.(round(primary.str.g.moderator.train.exercise.contrast[[2,3]],5))),
                                                                " [95% CI: ",bquote(.(round(primary.str.g.moderator.train.exercise.contrast[[2,6]],5))),
                                                                ", ",bquote(.(round(primary.str.g.moderator.train.exercise.contrast[[2,7]],5))),"]\n",
                                                             'Interaction Contrast (Multi-joint - Single-joint) = ',
                                                             bquote(.(round(primary.str.g.moderator.train.exercise.contrast[[3,3]],5))),
                                                             " [95% CI: ",bquote(.(round(primary.str.g.moderator.train.exercise.contrast[[3,6]],5))),
                                                             ", ",bquote(.(round(primary.str.g.moderator.train.exercise.contrast[[3,7]],5))),"]")


primary.str.g.moderator.train.exercise.plot=ggplot(data = primary.data.str%>%
                                                        mutate(train.exercise=ifelse(train.exercise=="Both","Both",
                                                                                     ifelse(train.exercise=="Multi-joint","Multi-joint Exercises","Single-joint Exercises"))),
                                                      aes(x=avg.rir,
                                                          y=yi))+THEME+theme(legend.position = "right",
                                                                             plot.caption = element_text(face = "italic"))+
  facet_wrap(~train.exercise)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Exercise Selection x Estimated RIR",
       caption = primary.str.g.moderator.train.exercise.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = primary.str.g.moderator.train.exercise.emms%>%
                mutate(train.exercise=ifelse(train.exercise=="Both","Both",
                                             ifelse(train.exercise=="Multi-joint","Multi-joint Exercises","Single-joint Exercises"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = primary.str.g.moderator.train.exercise.emms%>%
                mutate(train.exercise=ifelse(train.exercise=="Both","Both",
                                             ifelse(train.exercise=="Multi-joint","Multi-joint Exercises","Single-joint Exercises"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = primary.str.g.moderator.train.exercise.emms%>%
              mutate(train.exercise=ifelse(train.exercise=="Both","Both",
                                           ifelse(train.exercise=="Multi-joint","Multi-joint Exercises","Single-joint Exercises"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### g concurrent training
primary.str.g.moderator.formal.cardio.intervention=primary.str.g.best.model.prep%>%
  update(~. + log1p(avg.rir)*formal.cardio.intervention)%>%
  robust.rma.mv(cluster = primary.data.str$study)

primary.str.g.moderator.formal.cardio.intervention.emms=primary.str.g.moderator.formal.cardio.intervention%>%emmprep(at=list(avg.rir=0:23))%>%
  emmeans(~avg.rir+formal.cardio.intervention,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.str.g.moderator.formal.cardio.intervention)

primary.str.g.moderator.formal.cardio.intervention.contrast=emmprep(primary.str.g.moderator.formal.cardio.intervention,
                                                           at=list(avg.rir=primary.str.mean.rir))%>%
  emmeans(~avg.rir+formal.cardio.intervention,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()

primary.str.g.moderator.formal.cardio.intervention.contrast.label=paste0('Interaction Contrast (RT Only - Concurrent) = ',
                                                                bquote(.(round(primary.str.g.moderator.formal.cardio.intervention.contrast[[1,3]],5))),
                                                                " [95% CI: ",bquote(.(round(primary.str.g.moderator.formal.cardio.intervention.contrast[[1,6]],5))),
                                                                ", ",bquote(.(round(primary.str.g.moderator.formal.cardio.intervention.contrast[[1,7]],5))),"]")


primary.str.g.moderator.formal.cardio.intervention.plot=ggplot(data = primary.data.str%>%
                                                        mutate(formal.cardio.intervention=ifelse(formal.cardio.intervention=="No","Only Resistance Training","Concurrent Training")),
                                                      aes(x=avg.rir,
                                                          y=yi))+THEME+theme(legend.position = "right",
                                                                             plot.caption = element_text(face = "italic"))+
  facet_wrap(~formal.cardio.intervention)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Concurrent Training x Estimated RIR",
       caption = primary.str.g.moderator.formal.cardio.intervention.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = primary.str.g.moderator.formal.cardio.intervention.emms%>%
                mutate(formal.cardio.intervention=ifelse(formal.cardio.intervention=="No","Only Resistance Training","Concurrent Training")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = primary.str.g.moderator.formal.cardio.intervention.emms%>%
                mutate(formal.cardio.intervention=ifelse(formal.cardio.intervention=="No","Only Resistance Training","Concurrent Training")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = primary.str.g.moderator.formal.cardio.intervention.emms%>%
              mutate(formal.cardio.intervention=ifelse(formal.cardio.intervention=="No","Only Resistance Training","Concurrent Training")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### g progressive overload
primary.str.g.moderator.progression=primary.str.g.best.model.prep%>%
  update(~. + log1p(avg.rir)*progression)%>%
  robust.rma.mv(cluster = primary.data.str$study)

primary.str.g.moderator.progression.emms=primary.str.g.moderator.progression%>%emmprep(at=list(avg.rir=0:23))%>%
  emmeans(~avg.rir+progression,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.str.g.moderator.progression)

primary.str.g.moderator.progression.contrast=emmprep(primary.str.g.moderator.progression,
                                                           at=list(avg.rir=primary.str.mean.rir))%>%
  emmeans(~avg.rir+progression,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()

primary.str.g.moderator.progression.contrast.label=paste0('Interaction Contrast (No Overload - Progressive Overload) = ',
                                                                bquote(.(round(primary.str.g.moderator.progression.contrast[[1,3]],5))),
                                                                " [95% CI: ",bquote(.(round(primary.str.g.moderator.progression.contrast[[1,6]],5))),
                                                                ", ",bquote(.(round(primary.str.g.moderator.progression.contrast[[1,7]],5))),"]")


primary.str.g.moderator.progression.plot=ggplot(data = primary.data.str%>%
                                                        mutate(progression=ifelse(progression=="No","No Progressive Overload","Progressive Overload")),
                                                      aes(x=avg.rir,
                                                          y=yi))+THEME+theme(legend.position = "right",
                                                                             plot.caption = element_text(face = "italic"))+
  facet_wrap(~progression)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Progressive Overload x Estimated RIR",
       caption = primary.str.g.moderator.progression.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = primary.str.g.moderator.progression.emms%>%
                mutate(progression=ifelse(progression=="No","No Progressive Overload","Progressive Overload")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = primary.str.g.moderator.progression.emms%>%
                mutate(progression=ifelse(progression=="No","No Progressive Overload","Progressive Overload")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = primary.str.g.moderator.progression.emms%>%
              mutate(progression=ifelse(progression=="No","No Progressive Overload","Progressive Overload")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### g failure definition
primary.str.g.moderator.failure.definition=primary.str.g.best.model.prep%>%
  update(~. + log1p(avg.rir)*failure.definition)%>%
  robust.rma.mv(cluster = primary.data.str$study)

primary.str.g.moderator.failure.definition.emms=primary.str.g.moderator.failure.definition%>%emmprep(at=list(avg.rir=0:23))%>%
  emmeans(~avg.rir+failure.definition,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.str.g.moderator.failure.definition)

primary.str.g.moderator.failure.definition.contrast=emmprep(primary.str.g.moderator.failure.definition,
                                                           at=list(avg.rir=primary.str.mean.rir))%>%
  emmeans(~avg.rir+failure.definition,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()

primary.str.g.moderator.failure.definition.contrast.label=paste0('Interaction Contrast (No Definition - Definition Provided) = ',
                                                                bquote(.(round(primary.str.g.moderator.failure.definition.contrast[[1,3]],5))),
                                                                " [95% CI: ",bquote(.(round(primary.str.g.moderator.failure.definition.contrast[[1,6]],5))),
                                                                ", ",bquote(.(round(primary.str.g.moderator.failure.definition.contrast[[1,7]],5))),"]")


primary.str.g.moderator.failure.definition.plot=ggplot(data = primary.data.str%>%
                                                        mutate(failure.definition=ifelse(failure.definition=="No","No Failure Definition","Failure Definition Provided")),
                                                      aes(x=avg.rir,
                                                          y=yi))+THEME+theme(legend.position = "right",
                                                                             plot.caption = element_text(face = "italic"))+
  facet_wrap(~failure.definition)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Failure Definition x Estimated RIR",
       caption = primary.str.g.moderator.failure.definition.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = primary.str.g.moderator.failure.definition.emms%>%
                mutate(failure.definition=ifelse(failure.definition=="No","No Failure Definition","Failure Definition Provided")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = primary.str.g.moderator.failure.definition.emms%>%
                mutate(failure.definition=ifelse(failure.definition=="No","No Failure Definition","Failure Definition Provided")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = primary.str.g.moderator.failure.definition.emms%>%
              mutate(failure.definition=ifelse(failure.definition=="No","No Failure Definition","Failure Definition Provided")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### g set structure
primary.str.g.moderator.alternative.set.structure=primary.str.g.best.model.prep%>%
  update(~. + log1p(avg.rir)*alternative.set.structure)%>%
  robust.rma.mv(cluster = primary.data.str$study)

primary.str.g.moderator.alternative.set.structure.emms=primary.str.g.moderator.alternative.set.structure%>%emmprep(at=list(avg.rir=0:23))%>%
  emmeans(~avg.rir+alternative.set.structure,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.str.g.moderator.alternative.set.structure)

primary.str.g.moderator.alternative.set.structure.contrast=emmprep(primary.str.g.moderator.alternative.set.structure,
                                                           at=list(avg.rir=primary.str.mean.rir))%>%
  emmeans(~avg.rir+alternative.set.structure,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()

primary.str.g.moderator.alternative.set.structure.contrast.label=paste0('Interaction Contrast (Traditional - Alternative) = ',
                                                                bquote(.(round(primary.str.g.moderator.alternative.set.structure.contrast[[1,3]],5))),
                                                                " [95% CI: ",bquote(.(round(primary.str.g.moderator.alternative.set.structure.contrast[[1,6]],5))),
                                                                ", ",bquote(.(round(primary.str.g.moderator.alternative.set.structure.contrast[[1,7]],5))),"]")


primary.str.g.moderator.alternative.set.structure.plot=ggplot(data = primary.data.str%>%
                                                        mutate(alternative.set.structure=ifelse(alternative.set.structure=="No","Traditional Set Structure","Alternative Set Structure")),
                                                      aes(x=avg.rir,
                                                          y=yi))+THEME+theme(legend.position = "right",
                                                                             plot.caption = element_text(face = "italic"))+
  facet_wrap(~alternative.set.structure)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Set Structure x Estimated RIR",
       caption = primary.str.g.moderator.alternative.set.structure.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = primary.str.g.moderator.alternative.set.structure.emms%>%
                mutate(alternative.set.structure=ifelse(alternative.set.structure=="No","Traditional Set Structure","Alternative Set Structure")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = primary.str.g.moderator.alternative.set.structure.emms%>%
                mutate(alternative.set.structure=ifelse(alternative.set.structure=="No","Traditional Set Structure","Alternative Set Structure")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = primary.str.g.moderator.alternative.set.structure.emms%>%
              mutate(alternative.set.structure=ifelse(alternative.set.structure=="No","Traditional Set Structure","Alternative Set Structure")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### g intended velocity
primary.str.g.moderator.train.intent=primary.str.g.best.model.prep%>%
  update(~. + log1p(avg.rir)*train.intent)%>%
  robust.rma.mv(cluster = primary.data.str$study)

primary.str.g.moderator.train.intent.emms=primary.str.g.moderator.train.intent%>%emmprep(at=list(avg.rir=0:23))%>%
  emmeans(~avg.rir+train.intent,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.str.g.moderator.train.intent)

primary.str.g.moderator.train.intent.contrast=emmprep(primary.str.g.moderator.train.intent,
                                                           at=list(avg.rir=primary.str.mean.rir))%>%
  emmeans(~avg.rir+train.intent,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()

primary.str.g.moderator.train.intent.contrast.label=paste0('Interaction Contrast (Submaximal - Maximal) = ',
                                                                bquote(.(round(primary.str.g.moderator.train.intent.contrast[[1,3]],5))),
                                                                " [95% CI: ",bquote(.(round(primary.str.g.moderator.train.intent.contrast[[1,6]],5))),
                                                                ", ",bquote(.(round(primary.str.g.moderator.train.intent.contrast[[1,7]],5))),"]")


primary.str.g.moderator.train.intent.plot=ggplot(data = primary.data.str%>%
                                                        mutate(train.intent=ifelse(train.intent=="No","Submaximal Intended Concentric Velocity","Maximal Intended Concentric Velocity")),
                                                      aes(x=avg.rir,
                                                          y=yi))+THEME+theme(legend.position = "right",
                                                                             plot.caption = element_text(face = "italic"))+
  facet_wrap(~train.intent)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Intended Concentric Velocity x Estimated RIR",
       caption = primary.str.g.moderator.train.intent.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = primary.str.g.moderator.train.intent.emms%>%
                mutate(train.intent=ifelse(train.intent=="No","Submaximal Intended Concentric Velocity","Maximal Intended Concentric Velocity")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = primary.str.g.moderator.train.intent.emms%>%
                mutate(train.intent=ifelse(train.intent=="No","Submaximal Intended Concentric Velocity","Maximal Intended Concentric Velocity")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = primary.str.g.moderator.train.intent.emms%>%
              mutate(train.intent=ifelse(train.intent=="No","Submaximal Intended Concentric Velocity","Maximal Intended Concentric Velocity")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### g type of strength test
primary.str.g.moderator.isometric.isokinetic.isotonic=primary.str.g.best.model.prep%>%
  update(~. + log1p(avg.rir)*isometric.isokinetic.isotonic)%>%
  robust.rma.mv(cluster = primary.data.str$study)

primary.str.g.moderator.isometric.isokinetic.isotonic.emms=primary.str.g.moderator.isometric.isokinetic.isotonic%>%emmprep(at=list(avg.rir=0:23))%>%
  emmeans(~avg.rir+isometric.isokinetic.isotonic,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.str.g.moderator.isometric.isokinetic.isotonic)

primary.str.g.moderator.isometric.isokinetic.isotonic.contrast=emmprep(primary.str.g.moderator.isometric.isokinetic.isotonic,
                                                        at=list(avg.rir=primary.str.mean.rir))%>%
  emmeans(~avg.rir+isometric.isokinetic.isotonic,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()

primary.str.g.moderator.isometric.isokinetic.isotonic.contrast.label=paste0('Interaction Contrast (Isokinetic - Isometric) = ',
                                                             bquote(.(round(primary.str.g.moderator.isometric.isokinetic.isotonic.contrast[[1,3]],5))),
                                                             " [95% CI: ",bquote(.(round(primary.str.g.moderator.isometric.isokinetic.isotonic.contrast[[1,6]],5))),
                                                             ", ",bquote(.(round(primary.str.g.moderator.isometric.isokinetic.isotonic.contrast[[1,7]],5))),"]\n",
                                                             'Interaction Contrast (Isokinetic - Isotonic) = ',
                                                             bquote(.(round(primary.str.g.moderator.isometric.isokinetic.isotonic.contrast[[2,3]],5))),
                                                             " [95% CI: ",bquote(.(round(primary.str.g.moderator.isometric.isokinetic.isotonic.contrast[[2,6]],5))),
                                                             ", ",bquote(.(round(primary.str.g.moderator.isometric.isokinetic.isotonic.contrast[[2,7]],5))),"]\n",
                                                             'Interaction Contrast (Isometric - Isotonic) = ',
                                                             bquote(.(round(primary.str.g.moderator.isometric.isokinetic.isotonic.contrast[[3,3]],5))),
                                                             " [95% CI: ",bquote(.(round(primary.str.g.moderator.isometric.isokinetic.isotonic.contrast[[3,6]],5))),
                                                             ", ",bquote(.(round(primary.str.g.moderator.isometric.isokinetic.isotonic.contrast[[3,7]],5))),"]")


primary.str.g.moderator.isometric.isokinetic.isotonic.plot=ggplot(data = primary.data.str%>%
                                                     mutate(isometric.isokinetic.isotonic=ifelse(isometric.isokinetic.isotonic=="isokinetic","Isokinetic",
                                                                                  ifelse(isometric.isokinetic.isotonic=="isometric","Isometric","Isotonic"))),
                                                   aes(x=avg.rir,
                                                       y=yi))+THEME+theme(legend.position = "right",
                                                                          plot.caption = element_text(face = "italic"))+
  facet_wrap(~isometric.isokinetic.isotonic)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Type of Strength Test x Estimated RIR",
       caption = primary.str.g.moderator.isometric.isokinetic.isotonic.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = primary.str.g.moderator.isometric.isokinetic.isotonic.emms%>%
                mutate(isometric.isokinetic.isotonic=ifelse(isometric.isokinetic.isotonic=="isokinetic","Isokinetic",
                                                            ifelse(isometric.isokinetic.isotonic=="isometric","Isometric","Isotonic"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = primary.str.g.moderator.isometric.isokinetic.isotonic.emms%>%
                mutate(isometric.isokinetic.isotonic=ifelse(isometric.isokinetic.isotonic=="isokinetic","Isokinetic",
                                                            ifelse(isometric.isokinetic.isotonic=="isometric","Isometric","Isotonic"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = primary.str.g.moderator.isometric.isokinetic.isotonic.emms%>%
              mutate(isometric.isokinetic.isotonic=ifelse(isometric.isokinetic.isotonic=="isokinetic","Isokinetic",
                                                          ifelse(isometric.isokinetic.isotonic=="isometric","Isometric","Isotonic"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### g 1RM or non-1RM
primary.str.g.moderator.RM.max.submax=primary.str.g.best.model.prep%>%
  update(~. + log1p(avg.rir)*RM.max.submax)%>%
  robust.rma.mv(cluster = primary.data.str$study)

primary.str.g.moderator.RM.max.submax.emms=primary.str.g.moderator.RM.max.submax%>%emmprep(at=list(avg.rir=0:23))%>%
  emmeans(~avg.rir+RM.max.submax,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.str.g.moderator.RM.max.submax)

primary.str.g.moderator.RM.max.submax.contrast=emmprep(primary.str.g.moderator.RM.max.submax,
                                                           at=list(avg.rir=primary.str.mean.rir))%>%
  emmeans(~avg.rir+RM.max.submax,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()

primary.str.g.moderator.RM.max.submax.contrast.label=paste0('Interaction Contrast (1RM - Non 1RM) = ',
                                                                bquote(.(round(primary.str.g.moderator.RM.max.submax.contrast[[1,3]],5))),
                                                                " [95% CI: ",bquote(.(round(primary.str.g.moderator.RM.max.submax.contrast[[1,6]],5))),
                                                                ", ",bquote(.(round(primary.str.g.moderator.RM.max.submax.contrast[[1,7]],5))),"]")


primary.str.g.moderator.RM.max.submax.plot=ggplot(data = primary.data.str%>%
                                                        mutate(RM.max.submax=ifelse(RM.max.submax=="1RM","1RM","Non-1RM"))%>%
                                                    filter(RM.max.submax!="NA"),
                                                      aes(x=avg.rir,
                                                          y=yi))+THEME+theme(legend.position = "right",
                                                                             plot.caption = element_text(face = "italic"))+
  facet_wrap(~RM.max.submax)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Type of Isotonic Strength Test x Estimated RIR",
       caption = primary.str.g.moderator.RM.max.submax.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = primary.str.g.moderator.RM.max.submax.emms%>%
                mutate(RM.max.submax=ifelse(RM.max.submax=="1RM","1RM","Non-1RM")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = primary.str.g.moderator.RM.max.submax.emms%>%
                mutate(RM.max.submax=ifelse(RM.max.submax=="1RM","1RM","Non-1RM")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = primary.str.g.moderator.RM.max.submax.emms%>%
              mutate(RM.max.submax=ifelse(RM.max.submax=="1RM","1RM","Non-1RM")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### combine g plots
primary.str.g.moderator.doc=(primary.str.g.moderator.load.plot/plot_spacer()/primary.str.g.moderator.volume.method.plot/plot_spacer()/
  primary.str.g.moderator.train.status.plot/plot_spacer()/primary.str.g.moderator.weeks.plot/plot_spacer()/
  primary.str.g.moderator.adj.sets.plot/plot_spacer()/primary.str.g.moderator.frequency.per.muscle.plot/plot_spacer()/
  primary.str.g.moderator.age.plot/plot_spacer()/primary.str.g.moderator.sex.percent.male.plot/plot_spacer()/
  primary.str.g.moderator.within.between.design.plot/plot_spacer()/primary.str.g.moderator.upper.lower.other.plot/plot_spacer()/
  primary.str.g.moderator.train.exercise.plot/plot_spacer()/primary.str.g.moderator.formal.cardio.intervention.plot/plot_spacer()/
  primary.str.g.moderator.progression.plot/plot_spacer()/primary.str.g.moderator.failure.definition.plot/plot_spacer()/
  primary.str.g.moderator.alternative.set.structure.plot/plot_spacer()/primary.str.g.moderator.train.intent.plot/plot_spacer()/
  primary.str.g.moderator.isometric.isokinetic.isotonic.plot/plot_spacer()/primary.str.g.moderator.RM.max.submax.plot)&
  plot_annotation(title = "Moderator Analyses from Linear-log Multilevel Model for Maximal Strength (Standardized Mean Change)")&
  THEME

ggsave(primary.str.g.moderator.doc,device = "pdf",filename="primary.str.g.moderator.doc.pdf",width = 12,height = 156,limitsize = FALSE,
       path = "/Users/zacrobinson/Documents/For Publication PTF Meta/Plots")


#### rom load
primary.str.rom.moderator.load=primary.str.rom.best.model.prep%>%
  update(~. + avg.rir:load.set)%>%
  robust.rma.mv(cluster = primary.data.str$study)

primary.str.rom.moderator.load.emms=primary.str.rom.moderator.load%>%
  emmprep(at=list(avg.rir=0:23,
                  load.set=c(30,60,90)))%>%
  emmeans(~avg.rir+load.set,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.str.rom.moderator.load)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

primary.str.rom.moderator.load.contrast=emmprep(primary.str.rom.moderator.load,
                                              at=list(avg.rir=primary.str.mean.rir,
                                                      load.set=primary.str.mean.load.set))%>%
  emmeans(~avg.rir+load.set,
          weights = "prop")%>%
  contrast("consec",interaction=TRUE)%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

primary.str.rom.moderator.load.contrast.label=paste0('Interaction Contrast (Load) = ',
                                                   bquote(.(round(primary.str.rom.moderator.load.contrast[[3]],5))),
                                                   " [95% CI: ",bquote(.(round(primary.str.rom.moderator.load.contrast[[6]],5))),
                                                   ", ",bquote(.(round(primary.str.rom.moderator.load.contrast[[7]],5))),"]")

primary.str.rom.moderator.load.plot=ggplot(data = primary.data.str%>%
                                           mutate(load.set=ifelse(load.set<45,"30% of 1RM",
                                                                  ifelse(load.set>45&load.set<75,"60% of 1RM","90% of 1RM"))),
                                         aes(x=avg.rir,
                                             y=yi.rom.exp))+THEME+theme(plot.caption = element_text(face = "italic"))+
  facet_wrap(~load.set)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Load (% of 1RM) x Estimated RIR",
       caption = primary.str.rom.moderator.load.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = primary.str.rom.moderator.load.emms%>%
                mutate(load.set=ifelse(load.set==30,"30% of 1RM",
                                       ifelse(load.set==60,"60% of 1RM","90% of 1RM"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = primary.str.rom.moderator.load.emms%>%
                mutate(load.set=ifelse(load.set==30,"30% of 1RM",
                                       ifelse(load.set==60,"60% of 1RM","90% of 1RM"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = primary.str.rom.moderator.load.emms%>%
              mutate(load.set=ifelse(load.set==30,"30% of 1RM",
                                     ifelse(load.set==60,"60% of 1RM","90% of 1RM"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### rom volume equating method
primary.str.rom.moderator.volume.method=primary.str.rom.best.model.prep%>%
  update(~. + avg.rir:set.rep.equated)%>%
  robust.rma.mv(cluster = primary.data.str$study)

primary.str.rom.moderator.volume.method.emms=primary.str.rom.moderator.volume.method%>%emmprep(at=list(avg.rir=0:23))%>%
  emmeans(~avg.rir+set.rep.equated,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.str.rom.moderator.volume.method)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

primary.str.rom.moderator.volume.method.contrast=emmprep(primary.str.rom.moderator.volume.method,
                                                       at=list(avg.rir=primary.str.mean.rir))%>%
  emmeans(~avg.rir+set.rep.equated,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

primary.str.rom.moderator.volume.method.contrast.label=paste0('Interaction Contrast (Set - Rep) = ',
                                                            bquote(.(round(primary.str.rom.moderator.volume.method.contrast[[1,3]],5))),
                                                            " [95% CI: ",bquote(.(round(primary.str.rom.moderator.volume.method.contrast[[1,6]],5))),
                                                            ", ",bquote(.(round(primary.str.rom.moderator.volume.method.contrast[[1,7]],5))),"]\n",
                                                            'Interaction Contrast (Set - Both) = ',
                                                            bquote(.(round(primary.str.rom.moderator.volume.method.contrast[[2,3]],5))),
                                                            " [95% CI: ",bquote(.(round(primary.str.rom.moderator.volume.method.contrast[[2,6]],5))),
                                                            ", ",bquote(.(round(primary.str.rom.moderator.volume.method.contrast[[2,7]],5))),"]\n",
                                                            'Interaction Contrast (Rep - Both) = ',
                                                            bquote(.(round(primary.str.rom.moderator.volume.method.contrast[[3,3]],5))),
                                                            " [95% CI: ",bquote(.(round(primary.str.rom.moderator.volume.method.contrast[[3,6]],5))),
                                                            ", ",bquote(.(round(primary.str.rom.moderator.volume.method.contrast[[3,7]],5))),"]")


primary.str.rom.moderator.volume.method.plot=ggplot(data = primary.data.str%>%
                                                    mutate(set.rep.equated=ifelse(set.rep.equated=="set","Set Equated",
                                                                                  ifelse(set.rep.equated=="rep","Repetition Equated","Both Equated"))),
                                                  aes(x=avg.rir,
                                                      y=yi.rom.exp))+THEME+theme(legend.position = "right",
                                                                         plot.caption = element_text(face = "italic"))+
  facet_wrap(~set.rep.equated)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Method of Volume Equating x Estimated RIR",
       caption = primary.str.rom.moderator.volume.method.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = primary.str.rom.moderator.volume.method.emms%>%
                mutate(set.rep.equated=ifelse(set.rep.equated=="set","Set Equated",
                                              ifelse(set.rep.equated=="rep","Repetition Equated","Both Equated"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = primary.str.rom.moderator.volume.method.emms%>%
                mutate(set.rep.equated=ifelse(set.rep.equated=="set","Set Equated",
                                              ifelse(set.rep.equated=="rep","Repetition Equated","Both Equated"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = primary.str.rom.moderator.volume.method.emms%>%
              mutate(set.rep.equated=ifelse(set.rep.equated=="set","Set Equated",
                                            ifelse(set.rep.equated=="rep","Repetition Equated","Both Equated"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### rom training status
primary.str.rom.moderator.train.status=primary.str.rom.best.model.prep%>%
  update(~. + avg.rir:train.status)%>%
  robust.rma.mv(cluster = primary.data.str$study)

primary.str.rom.moderator.train.status.emms=primary.str.rom.moderator.train.status%>%emmprep(at=list(avg.rir=0:23))%>%
  emmeans(~avg.rir+train.status,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.str.rom.moderator.train.status)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

primary.str.rom.moderator.train.status.contrast=emmprep(primary.str.rom.moderator.train.status,
                                                      at=list(avg.rir=primary.str.mean.rir))%>%
  emmeans(~avg.rir+train.status,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

primary.str.rom.moderator.train.status.contrast.label=paste0('Interaction Contrast (Trained - Untrained) = ',
                                                           bquote(.(round(primary.str.rom.moderator.train.status.contrast[[1,3]],5))),
                                                           " [95% CI: ",bquote(.(round(primary.str.rom.moderator.train.status.contrast[[1,6]],5))),
                                                           ", ",bquote(.(round(primary.str.rom.moderator.train.status.contrast[[1,7]],5))),"]")


primary.str.rom.moderator.train.status.plot=ggplot(data = primary.data.str%>%
                                                   mutate(train.status=ifelse(train.status=="trained","Trained","Untrained")),
                                                 aes(x=avg.rir,
                                                     y=yi.rom.exp))+THEME+theme(legend.position = "right",
                                                                        plot.caption = element_text(face = "italic"))+
  facet_wrap(~train.status)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Training Status x Estimated RIR",
       caption = primary.str.rom.moderator.train.status.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = primary.str.rom.moderator.train.status.emms%>%
                mutate(train.status=ifelse(train.status=="trained","Trained","Untrained")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = primary.str.rom.moderator.train.status.emms%>%
                mutate(train.status=ifelse(train.status=="trained","Trained","Untrained")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = primary.str.rom.moderator.train.status.emms%>%
              mutate(train.status=ifelse(train.status=="trained","Trained","Untrained")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### rom weeks
primary.str.rom.moderator.weeks=primary.str.rom.best.model.prep%>%
  update(~. + avg.rir:weeks)%>%
  robust.rma.mv(cluster = primary.data.str$study)

primary.str.rom.moderator.weeks.emms=primary.str.rom.moderator.weeks%>%
  emmprep(at=list(avg.rir=0:23,
                  weeks=c(4,8,12)))%>%
  emmeans(~avg.rir+weeks,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.str.rom.moderator.weeks)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

primary.str.rom.moderator.weeks.contrast=emmprep(primary.str.rom.moderator.weeks,
                                               at=list(avg.rir=primary.str.mean.rir,
                                                       weeks=primary.str.mean.weeks))%>%
  emmeans(~avg.rir+weeks,
          weights = "prop")%>%
  contrast("consec",interaction=TRUE)%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

primary.str.rom.moderator.weeks.contrast.label=paste0('Interaction Contrast (Weeks) = ',
                                                    bquote(.(round(primary.str.rom.moderator.weeks.contrast[[3]],5))),
                                                    " [95% CI: ",bquote(.(round(primary.str.rom.moderator.weeks.contrast[[6]],5))),
                                                    ", ",bquote(.(round(primary.str.rom.moderator.weeks.contrast[[7]],5))),"]")

primary.str.rom.moderator.weeks.plot=ggplot(data = primary.data.str%>%
                                            mutate(weeks=ifelse(weeks<6,"4 Weeks",
                                                                ifelse(weeks>6&weeks<10,"8 Weeks","12 Weeks"))),
                                          aes(x=avg.rir,
                                              y=yi.rom.exp))+THEME+theme(plot.caption = element_text(face = "italic"))+
  facet_wrap(~fct_relevel(weeks,c("4 Weeks","8 Weeks","12 Weeks")))+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Intervention Duration (Weeks) x Estimated RIR",
       caption = primary.str.rom.moderator.weeks.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = primary.str.rom.moderator.weeks.emms%>%
                mutate(weeks=ifelse(weeks==4,"4 Weeks",
                                    ifelse(weeks==8,"8 Weeks","12 Weeks"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = primary.str.rom.moderator.weeks.emms%>%
                mutate(weeks=ifelse(weeks==4,"4 Weeks",
                                    ifelse(weeks==8,"8 Weeks","12 Weeks"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = primary.str.rom.moderator.weeks.emms%>%
              mutate(weeks=ifelse(weeks==4,"4 Weeks",
                                  ifelse(weeks==8,"8 Weeks","12 Weeks"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### rom adjusted sets 
primary.str.rom.moderator.adj.sets=primary.str.rom.best.model.prep%>%
  update(~. + avg.rir*adj.sets.week)%>%
  robust.rma.mv(cluster = primary.data.str$study)

primary.str.rom.moderator.adj.sets.emms=primary.str.rom.moderator.adj.sets%>%
  emmprep(at=list(avg.rir=0:23,
                  adj.sets.week=c(8,16,24)))%>%
  emmeans(~avg.rir+adj.sets.week,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.str.rom.moderator.adj.sets)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

primary.str.rom.moderator.adj.sets.contrast=emmprep(primary.str.rom.moderator.adj.sets,
                                                  at=list(avg.rir=primary.str.mean.rir,
                                                          adj.sets.week=primary.str.mean.adj.sets.week))%>%
  emmeans(~avg.rir+adj.sets.week,
          weights = "prop")%>%
  contrast("consec",interaction=TRUE)%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

primary.str.rom.moderator.adj.sets.contrast.label=paste0('Interaction Contrast (Adjusted Direct Sets) = ',
                                                       bquote(.(round(primary.str.rom.moderator.adj.sets.contrast[[3]],5))),
                                                       " [95% CI: ",bquote(.(round(primary.str.rom.moderator.adj.sets.contrast[[6]],5))),
                                                       ", ",bquote(.(round(primary.str.rom.moderator.adj.sets.contrast[[7]],5))),"]")

primary.str.rom.moderator.adj.sets.plot=ggplot(data = primary.data.str%>%
                                               mutate(adj.sets.week=ifelse(adj.sets.week<12,"8 Sets",
                                                                           ifelse(adj.sets.week>12&adj.sets.week<20,"16 Sets","24 Sets"))),
                                             aes(x=avg.rir,
                                                 y=yi.rom.exp))+THEME+theme(plot.caption = element_text(face = "italic"))+
  facet_wrap(~fct_relevel(adj.sets.week,c("8 Sets","16 Sets","24 Sets")))+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Adjusted Direct Sets Per Week Per Muscle Group x Estimated RIR",
       caption = primary.str.rom.moderator.adj.sets.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = primary.str.rom.moderator.adj.sets.emms%>%
                mutate(adj.sets.week=ifelse(adj.sets.week==8,"8 Sets",
                                            ifelse(adj.sets.week==16,"16 Sets","24 Sets"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = primary.str.rom.moderator.adj.sets.emms%>%
                mutate(adj.sets.week=ifelse(adj.sets.week==8,"8 Sets",
                                            ifelse(adj.sets.week==16,"16 Sets","24 Sets"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = primary.str.rom.moderator.adj.sets.emms%>%
              mutate(adj.sets.week=ifelse(adj.sets.week==8,"8 Sets",
                                          ifelse(adj.sets.week==16,"16 Sets","24 Sets"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### rom frequency
primary.str.rom.moderator.frequency.per.muscle=primary.str.rom.best.model.prep%>%
  update(~. + avg.rir*frequency.per.muscle)%>%
  robust.rma.mv(cluster = primary.data.str$study)

primary.str.rom.moderator.frequency.per.muscle.emms=primary.str.rom.moderator.frequency.per.muscle%>%
  emmprep(at=list(avg.rir=0:23,
                  frequency.per.muscle=c(1,2,3)))%>%
  emmeans(~avg.rir+frequency.per.muscle,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.str.rom.moderator.frequency.per.muscle)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

primary.str.rom.moderator.frequency.per.muscle.contrast=emmprep(primary.str.rom.moderator.frequency.per.muscle,
                                                              at=list(avg.rir=primary.str.mean.rir,
                                                                      frequency.per.muscle=primary.str.mean.frequency.per.muscle))%>%
  emmeans(~avg.rir+frequency.per.muscle,
          weights = "prop")%>%
  contrast("consec",interaction=TRUE)%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

primary.str.rom.moderator.frequency.per.muscle.contrast.label=paste0('Interaction Contrast (Direct Frequency) = ',
                                                                   bquote(.(round(primary.str.rom.moderator.frequency.per.muscle.contrast[[3]],5))),
                                                                   " [95% CI: ",bquote(.(round(primary.str.rom.moderator.frequency.per.muscle.contrast[[6]],5))),
                                                                   ", ",bquote(.(round(primary.str.rom.moderator.frequency.per.muscle.contrast[[7]],5))),"]")

primary.str.rom.moderator.frequency.per.muscle.plot=ggplot(data = primary.data.str%>%
                                                           mutate(frequency.per.muscle=ifelse(frequency.per.muscle<1.5,"1x Per Week",
                                                                                              ifelse(frequency.per.muscle>1.5&frequency.per.muscle<2.5,"2x Per Week","3x Per Week"))),
                                                         aes(x=avg.rir,
                                                             y=yi.rom.exp))+THEME+theme(plot.caption = element_text(face = "italic"))+
  facet_wrap(~frequency.per.muscle)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Direct Frequency Per Week Per Muscle Group x Estimated RIR",
       caption = primary.str.rom.moderator.frequency.per.muscle.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = primary.str.rom.moderator.frequency.per.muscle.emms%>%
                mutate(frequency.per.muscle=ifelse(frequency.per.muscle==1,"1x Per Week",
                                                   ifelse(frequency.per.muscle==2,"2x Per Week","3x Per Week"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = primary.str.rom.moderator.frequency.per.muscle.emms%>%
                mutate(frequency.per.muscle=ifelse(frequency.per.muscle==1,"1x Per Week",
                                                   ifelse(frequency.per.muscle==2,"2x Per Week","3x Per Week"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = primary.str.rom.moderator.frequency.per.muscle.emms%>%
              mutate(frequency.per.muscle=ifelse(frequency.per.muscle==1,"1x Per Week",
                                                 ifelse(frequency.per.muscle==2,"2x Per Week","3x Per Week"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### rom age
primary.str.rom.moderator.age=primary.str.rom.best.model.prep%>%
  update(~. + avg.rir*age)%>%
  robust.rma.mv(cluster = primary.data.str$study)

primary.str.rom.moderator.age.emms=primary.str.rom.moderator.age%>%
  emmprep(at=list(avg.rir=0:23,
                  age=c(20,40,60)))%>%
  emmeans(~avg.rir+age,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.str.rom.moderator.age)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

primary.str.rom.moderator.age.contrast=emmprep(primary.str.rom.moderator.age,
                                             at=list(avg.rir=primary.str.mean.rir,
                                                     age=primary.str.mean.age))%>%
  emmeans(~avg.rir+age,
          weights = "prop")%>%
  contrast("consec",interaction=TRUE)%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

primary.str.rom.moderator.age.contrast.label=paste0('Interaction Contrast (Age) = ',
                                                  bquote(.(round(primary.str.rom.moderator.age.contrast[[3]],5))),
                                                  " [95% CI: ",bquote(.(round(primary.str.rom.moderator.age.contrast[[6]],5))),
                                                  ", ",bquote(.(round(primary.str.rom.moderator.age.contrast[[7]],5))),"]")

primary.str.rom.moderator.age.plot=ggplot(data = primary.data.str%>%
                                          mutate(age=ifelse(age<30,"20 Years Old",
                                                            ifelse(age>30&age<50,"40 Years Old","60 Years Old"))),
                                        aes(x=avg.rir,
                                            y=yi.rom.exp))+THEME+theme(plot.caption = element_text(face = "italic"))+
  facet_wrap(~age)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Age x Estimated RIR",
       caption = primary.str.rom.moderator.age.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = primary.str.rom.moderator.age.emms%>%
                mutate(age=ifelse(age==20,"20 Years Old",
                                  ifelse(age==40,"40 Years Old","60 Years Old"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = primary.str.rom.moderator.age.emms%>%
                mutate(age=ifelse(age==20,"20 Years Old",
                                  ifelse(age==40,"40 Years Old","60 Years Old"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = primary.str.rom.moderator.age.emms%>%
              mutate(age=ifelse(age==20,"20 Years Old",
                                ifelse(age==40,"40 Years Old","60 Years Old"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### rom sex
primary.str.rom.moderator.sex.percent.male=primary.str.rom.best.model.prep%>%
  update(~. + avg.rir*sex.percent.male)%>%
  robust.rma.mv(cluster = primary.data.str$study)

primary.str.rom.moderator.sex.percent.male.emms=primary.str.rom.moderator.sex.percent.male%>%
  emmprep(at=list(avg.rir=0:23,
                  sex.percent.male=c(0,100)))%>%
  emmeans(~avg.rir+sex.percent.male,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.str.rom.moderator.sex.percent.male)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

primary.str.rom.moderator.sex.percent.male.contrast=emmprep(primary.str.rom.moderator.sex.percent.male,
                                                          at=list(avg.rir=primary.str.mean.rir,
                                                                  sex.percent.male=primary.str.mean.sex))%>%
  emmeans(~avg.rir+sex.percent.male,
          weights = "prop")%>%
  contrast("consec",interaction=TRUE)%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

primary.str.rom.moderator.sex.percent.male.contrast.label=paste0('Interaction Contrast (Sex) = ',
                                                               bquote(.(round(primary.str.rom.moderator.sex.percent.male.contrast[[3]],5))),
                                                               " [95% CI: ",bquote(.(round(primary.str.rom.moderator.sex.percent.male.contrast[[6]],5))),
                                                               ", ",bquote(.(round(primary.str.rom.moderator.sex.percent.male.contrast[[7]],5))),"]")

primary.str.rom.moderator.sex.percent.male.plot=ggplot(data = primary.data.str%>%
                                                       mutate(sex.percent.male=ifelse(sex.percent.male<50,"Female","Male")),
                                                     aes(x=avg.rir,
                                                         y=yi.rom.exp))+THEME+theme(plot.caption = element_text(face = "italic"))+
  facet_wrap(~sex.percent.male)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Sex (% Male) x Estimated RIR",
       caption = primary.str.rom.moderator.sex.percent.male.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = primary.str.rom.moderator.sex.percent.male.emms%>%
                mutate(sex.percent.male=ifelse(sex.percent.male==0,"Female","Male")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = primary.str.rom.moderator.sex.percent.male.emms%>%
                mutate(sex.percent.male=ifelse(sex.percent.male==0,"Female","Male")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = primary.str.rom.moderator.sex.percent.male.emms%>%
              mutate(sex.percent.male=ifelse(sex.percent.male==0,"Female","Male")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### rom study design
primary.str.rom.moderator.within.between.design=primary.str.rom.best.model.prep%>%
  update(~. + avg.rir*within.between.design)%>%
  robust.rma.mv(cluster = primary.data.str$study)

primary.str.rom.moderator.within.between.design.emms=primary.str.rom.moderator.within.between.design%>%emmprep(at=list(avg.rir=0:23))%>%
  emmeans(~avg.rir+within.between.design,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.str.rom.moderator.within.between.design)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

primary.str.rom.moderator.within.between.design.contrast=emmprep(primary.str.rom.moderator.within.between.design,
                                                               at=list(avg.rir=primary.str.mean.rir))%>%
  emmeans(~avg.rir+within.between.design,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

primary.str.rom.moderator.within.between.design.contrast.label=paste0('Interaction Contrast (Between - Within) = ',
                                                                    bquote(.(round(primary.str.rom.moderator.within.between.design.contrast[[1,3]],5))),
                                                                    " [95% CI: ",bquote(.(round(primary.str.rom.moderator.within.between.design.contrast[[1,6]],5))),
                                                                    ", ",bquote(.(round(primary.str.rom.moderator.within.between.design.contrast[[1,7]],5))),"]")


primary.str.rom.moderator.within.between.design.plot=ggplot(data = primary.data.str%>%
                                                            mutate(within.between.design=ifelse(within.between.design=="Between","Between Subjects","Within Subjects")),
                                                          aes(x=avg.rir,
                                                              y=yi.rom.exp))+THEME+theme(legend.position = "right",
                                                                                 plot.caption = element_text(face = "italic"))+
  facet_wrap(~within.between.design)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Study Design x Estimated RIR",
       caption = primary.str.rom.moderator.within.between.design.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = primary.str.rom.moderator.within.between.design.emms%>%
                mutate(within.between.design=ifelse(within.between.design=="Between","Between Subjects","Within Subjects")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = primary.str.rom.moderator.within.between.design.emms%>%
                mutate(within.between.design=ifelse(within.between.design=="Between","Between Subjects","Within Subjects")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = primary.str.rom.moderator.within.between.design.emms%>%
              mutate(within.between.design=ifelse(within.between.design=="Between","Between Subjects","Within Subjects")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### rom upper / lower outcome
primary.str.rom.moderator.upper.lower.other=primary.str.rom.best.model.prep%>%
  update(~. + avg.rir*upper.lower.other)%>%
  robust.rma.mv(cluster = primary.data.str$study)

primary.str.rom.moderator.upper.lower.other.emms=primary.str.rom.moderator.upper.lower.other%>%emmprep(at=list(avg.rir=0:23))%>%
  emmeans(~avg.rir+upper.lower.other,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.str.rom.moderator.upper.lower.other)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

primary.str.rom.moderator.upper.lower.other.contrast=emmprep(primary.str.rom.moderator.upper.lower.other,
                                                           at=list(avg.rir=primary.str.mean.rir))%>%
  emmeans(~avg.rir+upper.lower.other,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

primary.str.rom.moderator.upper.lower.other.contrast.label=paste0('Interaction Contrast (Lower - Upper) = ',
                                                                bquote(.(round(primary.str.rom.moderator.upper.lower.other.contrast[[1,3]],5))),
                                                                " [95% CI: ",bquote(.(round(primary.str.rom.moderator.upper.lower.other.contrast[[1,6]],5))),
                                                                ", ",bquote(.(round(primary.str.rom.moderator.upper.lower.other.contrast[[1,7]],5))),"]")


primary.str.rom.moderator.upper.lower.other.plot=ggplot(data = primary.data.str%>%
                                                        mutate(upper.lower.other=ifelse(upper.lower.other=="Lower","Lower Body Outcomes","Upper Body Outcomes")),
                                                      aes(x=avg.rir,
                                                          y=yi.rom.exp))+THEME+theme(legend.position = "right",
                                                                             plot.caption = element_text(face = "italic"))+
  facet_wrap(~upper.lower.other)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Upper/Lower Body Outcome x Estimated RIR",
       caption = primary.str.rom.moderator.upper.lower.other.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = primary.str.rom.moderator.upper.lower.other.emms%>%
                mutate(upper.lower.other=ifelse(upper.lower.other=="Lower","Lower Body Outcomes","Upper Body Outcomes")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = primary.str.rom.moderator.upper.lower.other.emms%>%
                mutate(upper.lower.other=ifelse(upper.lower.other=="Lower","Lower Body Outcomes","Upper Body Outcomes")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = primary.str.rom.moderator.upper.lower.other.emms%>%
              mutate(upper.lower.other=ifelse(upper.lower.other=="Lower","Lower Body Outcomes","Upper Body Outcomes")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### rom exercise selection
primary.str.rom.moderator.train.exercise=primary.str.rom.best.model.prep%>%
  update(~. + avg.rir*train.exercise)%>%
  robust.rma.mv(cluster = primary.data.str$study)

primary.str.rom.moderator.train.exercise.emms=primary.str.rom.moderator.train.exercise%>%emmprep(at=list(avg.rir=0:23))%>%
  emmeans(~avg.rir+train.exercise,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.str.rom.moderator.train.exercise)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

primary.str.rom.moderator.train.exercise.contrast=emmprep(primary.str.rom.moderator.train.exercise,
                                                        at=list(avg.rir=primary.str.mean.rir))%>%
  emmeans(~avg.rir+train.exercise,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

primary.str.rom.moderator.train.exercise.contrast.label=paste0('Interaction Contrast (Both - Multi-joint) = ',
                                                             bquote(.(round(primary.str.rom.moderator.train.exercise.contrast[[1,3]],5))),
                                                             " [95% CI: ",bquote(.(round(primary.str.rom.moderator.train.exercise.contrast[[1,6]],5))),
                                                             ", ",bquote(.(round(primary.str.rom.moderator.train.exercise.contrast[[1,7]],5))),"]\n",
                                                             'Interaction Contrast (Both - Single-joint) = ',
                                                             bquote(.(round(primary.str.rom.moderator.train.exercise.contrast[[2,3]],5))),
                                                             " [95% CI: ",bquote(.(round(primary.str.rom.moderator.train.exercise.contrast[[2,6]],5))),
                                                             ", ",bquote(.(round(primary.str.rom.moderator.train.exercise.contrast[[2,7]],5))),"]\n",
                                                             'Interaction Contrast (Multi-joint - Single-joint) = ',
                                                             bquote(.(round(primary.str.rom.moderator.train.exercise.contrast[[3,3]],5))),
                                                             " [95% CI: ",bquote(.(round(primary.str.rom.moderator.train.exercise.contrast[[3,6]],5))),
                                                             ", ",bquote(.(round(primary.str.rom.moderator.train.exercise.contrast[[3,7]],5))),"]")


primary.str.rom.moderator.train.exercise.plot=ggplot(data = primary.data.str%>%
                                                     mutate(train.exercise=ifelse(train.exercise=="Both","Both",
                                                                                  ifelse(train.exercise=="Multi-joint","Multi-joint Exercises","Single-joint Exercises"))),
                                                   aes(x=avg.rir,
                                                       y=yi.rom.exp))+THEME+theme(legend.position = "right",
                                                                          plot.caption = element_text(face = "italic"))+
  facet_wrap(~train.exercise)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Exercise Selection x Estimated RIR",
       caption = primary.str.rom.moderator.train.exercise.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = primary.str.rom.moderator.train.exercise.emms%>%
                mutate(train.exercise=ifelse(train.exercise=="Both","Both",
                                             ifelse(train.exercise=="Multi-joint","Multi-joint Exercises","Single-joint Exercises"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = primary.str.rom.moderator.train.exercise.emms%>%
                mutate(train.exercise=ifelse(train.exercise=="Both","Both",
                                             ifelse(train.exercise=="Multi-joint","Multi-joint Exercises","Single-joint Exercises"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = primary.str.rom.moderator.train.exercise.emms%>%
              mutate(train.exercise=ifelse(train.exercise=="Both","Both",
                                           ifelse(train.exercise=="Multi-joint","Multi-joint Exercises","Single-joint Exercises"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### rom concurrent training
primary.str.rom.moderator.formal.cardio.intervention=primary.str.rom.best.model.prep%>%
  update(~. + avg.rir*formal.cardio.intervention)%>%
  robust.rma.mv(cluster = primary.data.str$study)

primary.str.rom.moderator.formal.cardio.intervention.emms=primary.str.rom.moderator.formal.cardio.intervention%>%emmprep(at=list(avg.rir=0:23))%>%
  emmeans(~avg.rir+formal.cardio.intervention,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.str.rom.moderator.formal.cardio.intervention)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

primary.str.rom.moderator.formal.cardio.intervention.contrast=emmprep(primary.str.rom.moderator.formal.cardio.intervention,
                                                                    at=list(avg.rir=primary.str.mean.rir))%>%
  emmeans(~avg.rir+formal.cardio.intervention,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

primary.str.rom.moderator.formal.cardio.intervention.contrast.label=paste0('Interaction Contrast (RT Only - Concurrent) = ',
                                                                         bquote(.(round(primary.str.rom.moderator.formal.cardio.intervention.contrast[[1,3]],5))),
                                                                         " [95% CI: ",bquote(.(round(primary.str.rom.moderator.formal.cardio.intervention.contrast[[1,6]],5))),
                                                                         ", ",bquote(.(round(primary.str.rom.moderator.formal.cardio.intervention.contrast[[1,7]],5))),"]")


primary.str.rom.moderator.formal.cardio.intervention.plot=ggplot(data = primary.data.str%>%
                                                                 mutate(formal.cardio.intervention=ifelse(formal.cardio.intervention=="No","Only Resistance Training","Concurrent Training")),
                                                               aes(x=avg.rir,
                                                                   y=yi.rom.exp))+THEME+theme(legend.position = "right",
                                                                                      plot.caption = element_text(face = "italic"))+
  facet_wrap(~formal.cardio.intervention)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Concurrent Training x Estimated RIR",
       caption = primary.str.rom.moderator.formal.cardio.intervention.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = primary.str.rom.moderator.formal.cardio.intervention.emms%>%
                mutate(formal.cardio.intervention=ifelse(formal.cardio.intervention=="No","Only Resistance Training","Concurrent Training")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = primary.str.rom.moderator.formal.cardio.intervention.emms%>%
                mutate(formal.cardio.intervention=ifelse(formal.cardio.intervention=="No","Only Resistance Training","Concurrent Training")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = primary.str.rom.moderator.formal.cardio.intervention.emms%>%
              mutate(formal.cardio.intervention=ifelse(formal.cardio.intervention=="No","Only Resistance Training","Concurrent Training")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### rom progressive overload
primary.str.rom.moderator.progression=primary.str.rom.best.model.prep%>%
  update(~. + avg.rir*progression)%>%
  robust.rma.mv(cluster = primary.data.str$study)

primary.str.rom.moderator.progression.emms=primary.str.rom.moderator.progression%>%emmprep(at=list(avg.rir=0:23))%>%
  emmeans(~avg.rir+progression,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.str.rom.moderator.progression)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

primary.str.rom.moderator.progression.contrast=emmprep(primary.str.rom.moderator.progression,
                                                     at=list(avg.rir=primary.str.mean.rir))%>%
  emmeans(~avg.rir+progression,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

primary.str.rom.moderator.progression.contrast.label=paste0('Interaction Contrast (No Overload - Progressive Overload) = ',
                                                          bquote(.(round(primary.str.rom.moderator.progression.contrast[[1,3]],5))),
                                                          " [95% CI: ",bquote(.(round(primary.str.rom.moderator.progression.contrast[[1,6]],5))),
                                                          ", ",bquote(.(round(primary.str.rom.moderator.progression.contrast[[1,7]],5))),"]")


primary.str.rom.moderator.progression.plot=ggplot(data = primary.data.str%>%
                                                  mutate(progression=ifelse(progression=="No","No Progressive Overload","Progressive Overload")),
                                                aes(x=avg.rir,
                                                    y=yi.rom.exp))+THEME+theme(legend.position = "right",
                                                                       plot.caption = element_text(face = "italic"))+
  facet_wrap(~progression)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Progressive Overload x Estimated RIR",
       caption = primary.str.rom.moderator.progression.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = primary.str.rom.moderator.progression.emms%>%
                mutate(progression=ifelse(progression=="No","No Progressive Overload","Progressive Overload")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = primary.str.rom.moderator.progression.emms%>%
                mutate(progression=ifelse(progression=="No","No Progressive Overload","Progressive Overload")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = primary.str.rom.moderator.progression.emms%>%
              mutate(progression=ifelse(progression=="No","No Progressive Overload","Progressive Overload")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### rom failure definition
primary.str.rom.moderator.failure.definition=primary.str.rom.best.model.prep%>%
  update(~. + avg.rir*failure.definition)%>%
  robust.rma.mv(cluster = primary.data.str$study)

primary.str.rom.moderator.failure.definition.emms=primary.str.rom.moderator.failure.definition%>%emmprep(at=list(avg.rir=0:23))%>%
  emmeans(~avg.rir+failure.definition,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.str.rom.moderator.failure.definition)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

primary.str.rom.moderator.failure.definition.contrast=emmprep(primary.str.rom.moderator.failure.definition,
                                                            at=list(avg.rir=primary.str.mean.rir))%>%
  emmeans(~avg.rir+failure.definition,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

primary.str.rom.moderator.failure.definition.contrast.label=paste0('Interaction Contrast (No Definition - Definition Provided) = ',
                                                                 bquote(.(round(primary.str.rom.moderator.failure.definition.contrast[[1,3]],5))),
                                                                 " [95% CI: ",bquote(.(round(primary.str.rom.moderator.failure.definition.contrast[[1,6]],5))),
                                                                 ", ",bquote(.(round(primary.str.rom.moderator.failure.definition.contrast[[1,7]],5))),"]")


primary.str.rom.moderator.failure.definition.plot=ggplot(data = primary.data.str%>%
                                                         mutate(failure.definition=ifelse(failure.definition=="No","No Failure Definition","Failure Definition Provided")),
                                                       aes(x=avg.rir,
                                                           y=yi.rom.exp))+THEME+theme(legend.position = "right",
                                                                              plot.caption = element_text(face = "italic"))+
  facet_wrap(~failure.definition)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Failure Definition x Estimated RIR",
       caption = primary.str.rom.moderator.failure.definition.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = primary.str.rom.moderator.failure.definition.emms%>%
                mutate(failure.definition=ifelse(failure.definition=="No","No Failure Definition","Failure Definition Provided")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = primary.str.rom.moderator.failure.definition.emms%>%
                mutate(failure.definition=ifelse(failure.definition=="No","No Failure Definition","Failure Definition Provided")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = primary.str.rom.moderator.failure.definition.emms%>%
              mutate(failure.definition=ifelse(failure.definition=="No","No Failure Definition","Failure Definition Provided")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### rom set structure
primary.str.rom.moderator.alternative.set.structure=primary.str.rom.best.model.prep%>%
  update(~. + avg.rir*alternative.set.structure)%>%
  robust.rma.mv(cluster = primary.data.str$study)

primary.str.rom.moderator.alternative.set.structure.emms=primary.str.rom.moderator.alternative.set.structure%>%emmprep(at=list(avg.rir=0:23))%>%
  emmeans(~avg.rir+alternative.set.structure,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.str.rom.moderator.alternative.set.structure)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

primary.str.rom.moderator.alternative.set.structure.contrast=emmprep(primary.str.rom.moderator.alternative.set.structure,
                                                                   at=list(avg.rir=primary.str.mean.rir))%>%
  emmeans(~avg.rir+alternative.set.structure,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

primary.str.rom.moderator.alternative.set.structure.contrast.label=paste0('Interaction Contrast (Traditional - Alternative) = ',
                                                                        bquote(.(round(primary.str.rom.moderator.alternative.set.structure.contrast[[1,3]],5))),
                                                                        " [95% CI: ",bquote(.(round(primary.str.rom.moderator.alternative.set.structure.contrast[[1,6]],5))),
                                                                        ", ",bquote(.(round(primary.str.rom.moderator.alternative.set.structure.contrast[[1,7]],5))),"]")


primary.str.rom.moderator.alternative.set.structure.plot=ggplot(data = primary.data.str%>%
                                                                mutate(alternative.set.structure=ifelse(alternative.set.structure=="No","Traditional Set Structure","Alternative Set Structure")),
                                                              aes(x=avg.rir,
                                                                  y=yi.rom.exp))+THEME+theme(legend.position = "right",
                                                                                     plot.caption = element_text(face = "italic"))+
  facet_wrap(~alternative.set.structure)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Set Structure x Estimated RIR",
       caption = primary.str.rom.moderator.alternative.set.structure.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = primary.str.rom.moderator.alternative.set.structure.emms%>%
                mutate(alternative.set.structure=ifelse(alternative.set.structure=="No","Traditional Set Structure","Alternative Set Structure")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = primary.str.rom.moderator.alternative.set.structure.emms%>%
                mutate(alternative.set.structure=ifelse(alternative.set.structure=="No","Traditional Set Structure","Alternative Set Structure")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = primary.str.rom.moderator.alternative.set.structure.emms%>%
              mutate(alternative.set.structure=ifelse(alternative.set.structure=="No","Traditional Set Structure","Alternative Set Structure")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### rom intended velocity
primary.str.rom.moderator.train.intent=primary.str.rom.best.model.prep%>%
  update(~. + avg.rir*train.intent)%>%
  robust.rma.mv(cluster = primary.data.str$study)

primary.str.rom.moderator.train.intent.emms=primary.str.rom.moderator.train.intent%>%emmprep(at=list(avg.rir=0:23))%>%
  emmeans(~avg.rir+train.intent,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.str.rom.moderator.train.intent)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

primary.str.rom.moderator.train.intent.contrast=emmprep(primary.str.rom.moderator.train.intent,
                                                      at=list(avg.rir=primary.str.mean.rir))%>%
  emmeans(~avg.rir+train.intent,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

primary.str.rom.moderator.train.intent.contrast.label=paste0('Interaction Contrast (Submaximal - Maximal) = ',
                                                           bquote(.(round(primary.str.rom.moderator.train.intent.contrast[[1,3]],5))),
                                                           " [95% CI: ",bquote(.(round(primary.str.rom.moderator.train.intent.contrast[[1,6]],5))),
                                                           ", ",bquote(.(round(primary.str.rom.moderator.train.intent.contrast[[1,7]],5))),"]")


primary.str.rom.moderator.train.intent.plot=ggplot(data = primary.data.str%>%
                                                   mutate(train.intent=ifelse(train.intent=="No","Submaximal Intended Concentric Velocity","Maximal Intended Concentric Velocity")),
                                                 aes(x=avg.rir,
                                                     y=yi.rom.exp))+THEME+theme(legend.position = "right",
                                                                        plot.caption = element_text(face = "italic"))+
  facet_wrap(~train.intent)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Intended Concentric Velocity x Estimated RIR",
       caption = primary.str.rom.moderator.train.intent.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = primary.str.rom.moderator.train.intent.emms%>%
                mutate(train.intent=ifelse(train.intent=="No","Submaximal Intended Concentric Velocity","Maximal Intended Concentric Velocity")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = primary.str.rom.moderator.train.intent.emms%>%
                mutate(train.intent=ifelse(train.intent=="No","Submaximal Intended Concentric Velocity","Maximal Intended Concentric Velocity")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = primary.str.rom.moderator.train.intent.emms%>%
              mutate(train.intent=ifelse(train.intent=="No","Submaximal Intended Concentric Velocity","Maximal Intended Concentric Velocity")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### rom type of strength test
primary.str.rom.moderator.isometric.isokinetic.isotonic=primary.str.rom.best.model.prep%>%
  update(~. + avg.rir*isometric.isokinetic.isotonic)%>%
  robust.rma.mv(cluster = primary.data.str$study)

primary.str.rom.moderator.isometric.isokinetic.isotonic.emms=primary.str.rom.moderator.isometric.isokinetic.isotonic%>%emmprep(at=list(avg.rir=0:23))%>%
  emmeans(~avg.rir+isometric.isokinetic.isotonic,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.str.rom.moderator.isometric.isokinetic.isotonic)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

primary.str.rom.moderator.isometric.isokinetic.isotonic.contrast=emmprep(primary.str.rom.moderator.isometric.isokinetic.isotonic,
                                                                       at=list(avg.rir=primary.str.mean.rir))%>%
  emmeans(~avg.rir+isometric.isokinetic.isotonic,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

primary.str.rom.moderator.isometric.isokinetic.isotonic.contrast.label=paste0('Interaction Contrast (Isokinetic - Isometric) = ',
                                                                            bquote(.(round(primary.str.rom.moderator.isometric.isokinetic.isotonic.contrast[[1,3]],5))),
                                                                            " [95% CI: ",bquote(.(round(primary.str.rom.moderator.isometric.isokinetic.isotonic.contrast[[1,6]],5))),
                                                                            ", ",bquote(.(round(primary.str.rom.moderator.isometric.isokinetic.isotonic.contrast[[1,7]],5))),"]\n",
                                                                            'Interaction Contrast (Isokinetic - Isotonic) = ',
                                                                            bquote(.(round(primary.str.rom.moderator.isometric.isokinetic.isotonic.contrast[[2,3]],5))),
                                                                            " [95% CI: ",bquote(.(round(primary.str.rom.moderator.isometric.isokinetic.isotonic.contrast[[2,6]],5))),
                                                                            ", ",bquote(.(round(primary.str.rom.moderator.isometric.isokinetic.isotonic.contrast[[2,7]],5))),"]\n",
                                                                            'Interaction Contrast (Isometric - Isotonic) = ',
                                                                            bquote(.(round(primary.str.rom.moderator.isometric.isokinetic.isotonic.contrast[[3,3]],5))),
                                                                            " [95% CI: ",bquote(.(round(primary.str.rom.moderator.isometric.isokinetic.isotonic.contrast[[3,6]],5))),
                                                                            ", ",bquote(.(round(primary.str.rom.moderator.isometric.isokinetic.isotonic.contrast[[3,7]],5))),"]")


primary.str.rom.moderator.isometric.isokinetic.isotonic.plot=ggplot(data = primary.data.str%>%
                                                                    mutate(isometric.isokinetic.isotonic=ifelse(isometric.isokinetic.isotonic=="isokinetic","Isokinetic",
                                                                                                                ifelse(isometric.isokinetic.isotonic=="isometric","Isometric","Isotonic"))),
                                                                  aes(x=avg.rir,
                                                                      y=yi.rom.exp))+THEME+theme(legend.position = "right",
                                                                                         plot.caption = element_text(face = "italic"))+
  facet_wrap(~isometric.isokinetic.isotonic)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Type of Strength Test x Estimated RIR",
       caption = primary.str.rom.moderator.isometric.isokinetic.isotonic.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = primary.str.rom.moderator.isometric.isokinetic.isotonic.emms%>%
                mutate(isometric.isokinetic.isotonic=ifelse(isometric.isokinetic.isotonic=="isokinetic","Isokinetic",
                                                            ifelse(isometric.isokinetic.isotonic=="isometric","Isometric","Isotonic"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = primary.str.rom.moderator.isometric.isokinetic.isotonic.emms%>%
                mutate(isometric.isokinetic.isotonic=ifelse(isometric.isokinetic.isotonic=="isokinetic","Isokinetic",
                                                            ifelse(isometric.isokinetic.isotonic=="isometric","Isometric","Isotonic"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = primary.str.rom.moderator.isometric.isokinetic.isotonic.emms%>%
              mutate(isometric.isokinetic.isotonic=ifelse(isometric.isokinetic.isotonic=="isokinetic","Isokinetic",
                                                          ifelse(isometric.isokinetic.isotonic=="isometric","Isometric","Isotonic"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### rom 1RM or non-1RM
primary.str.rom.moderator.RM.max.submax=primary.str.rom.best.model.prep%>%
  update(~. + avg.rir*RM.max.submax)%>%
  robust.rma.mv(cluster = primary.data.str$study)

primary.str.rom.moderator.RM.max.submax.emms=primary.str.rom.moderator.RM.max.submax%>%emmprep(at=list(avg.rir=0:23))%>%
  emmeans(~avg.rir+RM.max.submax,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.str.rom.moderator.RM.max.submax)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

primary.str.rom.moderator.RM.max.submax.contrast=emmprep(primary.str.rom.moderator.RM.max.submax,
                                                       at=list(avg.rir=primary.str.mean.rir))%>%
  emmeans(~avg.rir+RM.max.submax,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

primary.str.rom.moderator.RM.max.submax.contrast.label=paste0('Interaction Contrast (1RM - Non 1RM) = ',
                                                            bquote(.(round(primary.str.rom.moderator.RM.max.submax.contrast[[1,3]],5))),
                                                            " [95% CI: ",bquote(.(round(primary.str.rom.moderator.RM.max.submax.contrast[[1,6]],5))),
                                                            ", ",bquote(.(round(primary.str.rom.moderator.RM.max.submax.contrast[[1,7]],5))),"]")


primary.str.rom.moderator.RM.max.submax.plot=ggplot(data = primary.data.str%>%
                                                    mutate(RM.max.submax=ifelse(RM.max.submax=="1RM","1RM","Non-1RM"))%>%
                                                    filter(RM.max.submax!="NA"),
                                                  aes(x=avg.rir,
                                                      y=yi.rom.exp))+THEME+theme(legend.position = "right",
                                                                         plot.caption = element_text(face = "italic"))+
  facet_wrap(~RM.max.submax)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Type of Isotonic Strength Test x Estimated RIR",
       caption = primary.str.rom.moderator.RM.max.submax.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = primary.str.rom.moderator.RM.max.submax.emms%>%
                mutate(RM.max.submax=ifelse(RM.max.submax=="1RM","1RM","Non-1RM")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = primary.str.rom.moderator.RM.max.submax.emms%>%
                mutate(RM.max.submax=ifelse(RM.max.submax=="1RM","1RM","Non-1RM")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = primary.str.rom.moderator.RM.max.submax.emms%>%
              mutate(RM.max.submax=ifelse(RM.max.submax=="1RM","1RM","Non-1RM")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### combine rom plots
primary.str.rom.moderator.doc=primary.str.rom.moderator.load.plot/plot_spacer()/primary.str.rom.moderator.volume.method.plot/plot_spacer()/
  primary.str.rom.moderator.train.status.plot/plot_spacer()/primary.str.rom.moderator.weeks.plot/plot_spacer()/
  primary.str.rom.moderator.adj.sets.plot/plot_spacer()/primary.str.rom.moderator.frequency.per.muscle.plot/plot_spacer()/
  primary.str.rom.moderator.age.plot/plot_spacer()/primary.str.rom.moderator.sex.percent.male.plot/plot_spacer()/
  primary.str.rom.moderator.within.between.design.plot/plot_spacer()/primary.str.rom.moderator.upper.lower.other.plot/plot_spacer()/
  primary.str.rom.moderator.train.exercise.plot/plot_spacer()/primary.str.rom.moderator.formal.cardio.intervention.plot/plot_spacer()/
  primary.str.rom.moderator.progression.plot/plot_spacer()/primary.str.rom.moderator.failure.definition.plot/plot_spacer()/
  primary.str.rom.moderator.alternative.set.structure.plot/plot_spacer()/primary.str.rom.moderator.train.intent.plot/plot_spacer()/
  primary.str.rom.moderator.isometric.isokinetic.isotonic.plot/plot_spacer()/primary.str.rom.moderator.RM.max.submax.plot&
  plot_annotation(title = "Moderator Analyses from Linear Fit Multilevel Model for Maximal Strength (Response Ratio)")&
  THEME

ggsave(primary.str.rom.moderator.doc,device = "pdf",filename="primary.str.rom.moderator.doc.pdf",width = 12,height = 156,limitsize = FALSE,
       path = "/Users/zacrobinson/Documents/For Publication PTF Meta/Plots")

### hypertrophy
#### g load
primary.hyp.g.moderator.load=primary.hyp.g.best.model.prep%>%
  update(~. + avg.rir:load.set)%>%
  robust.rma.mv(cluster = primary.data.hyp$study)

primary.hyp.g.moderator.load.emms=primary.hyp.g.moderator.load%>%
  emmprep(at=list(avg.rir=0:23,
                  load.set=c(30,60,90)))%>%
  emmeans(~avg.rir+load.set,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.hyp.g.moderator.load)

primary.hyp.g.moderator.load.contrast=emmprep(primary.hyp.g.moderator.load,
                                              at=list(avg.rir=primary.hyp.mean.rir,
                                                      load.set=primary.hyp.mean.load.set))%>%
  emmeans(~avg.rir+load.set,
          weights = "prop")%>%
  contrast(interaction= c("consec","pairwise"))%>%
  confint()

primary.hyp.g.moderator.load.contrast.label=paste0('Interaction Contrast (Load) = ',
                                                   bquote(.(round(primary.hyp.g.moderator.load.contrast[[3]],5))),
                                                   " [95% CI: ",bquote(.(round(primary.hyp.g.moderator.load.contrast[[6]],5))),
                                                   ", ",bquote(.(round(primary.hyp.g.moderator.load.contrast[[7]],5))),"]")

primary.hyp.g.moderator.load.plot=ggplot(data = primary.data.hyp%>%
                                           mutate(load.set=ifelse(load.set<45,"30% of 1RM",
                                                                  ifelse(load.set>45&load.set<75,"60% of 1RM","90% of 1RM"))),
                                         aes(x=avg.rir,
                                             y=yi))+THEME+theme(plot.caption = element_text(face = "italic"))+
  facet_wrap(~load.set)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Load (% of 1RM) x Estimated RIR",
       caption = primary.hyp.g.moderator.load.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = primary.hyp.g.moderator.load.emms%>%
                mutate(load.set=ifelse(load.set==30,"30% of 1RM",
                                       ifelse(load.set==60,"60% of 1RM","90% of 1RM"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = primary.hyp.g.moderator.load.emms%>%
                mutate(load.set=ifelse(load.set==30,"30% of 1RM",
                                       ifelse(load.set==60,"60% of 1RM","90% of 1RM"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = primary.hyp.g.moderator.load.emms%>%
              mutate(load.set=ifelse(load.set==30,"30% of 1RM",
                                     ifelse(load.set==60,"60% of 1RM","90% of 1RM"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### g volume equating method
primary.hyp.g.moderator.volume.method=primary.hyp.g.best.model.prep%>%
  update(~. + avg.rir:set.rep.equated)%>%
  robust.rma.mv(cluster = primary.data.hyp$study)

primary.hyp.g.moderator.volume.method.emms=primary.hyp.g.moderator.volume.method%>%emmprep(at=list(avg.rir=0:23))%>%
  emmeans(~avg.rir+set.rep.equated,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.hyp.g.moderator.volume.method)

primary.hyp.g.moderator.volume.method.contrast=emmprep(primary.hyp.g.moderator.volume.method,
                                                       at=list(avg.rir=primary.hyp.mean.rir))%>%
  emmeans(~avg.rir+set.rep.equated,
          weights = "prop")%>%
  contrast(interaction= c("consec","pairwise"))%>%
  confint()

primary.hyp.g.moderator.volume.method.contrast.label=paste0('Interaction Contrast (Set - Rep) = ',
                                                            bquote(.(round(primary.hyp.g.moderator.volume.method.contrast[[1,3]],5))),
                                                            " [95% CI: ",bquote(.(round(primary.hyp.g.moderator.volume.method.contrast[[1,6]],5))),
                                                            ", ",bquote(.(round(primary.hyp.g.moderator.volume.method.contrast[[1,7]],5))),"]\n",
                                                            'Interaction Contrast (Set - Both) = ',
                                                            bquote(.(round(primary.hyp.g.moderator.volume.method.contrast[[2,3]],5))),
                                                            " [95% CI: ",bquote(.(round(primary.hyp.g.moderator.volume.method.contrast[[2,6]],5))),
                                                            ", ",bquote(.(round(primary.hyp.g.moderator.volume.method.contrast[[2,7]],5))),"]\n",
                                                            'Interaction Contrast (Rep - Both) = ',
                                                            bquote(.(round(primary.hyp.g.moderator.volume.method.contrast[[3,3]],5))),
                                                            " [95% CI: ",bquote(.(round(primary.hyp.g.moderator.volume.method.contrast[[3,6]],5))),
                                                            ", ",bquote(.(round(primary.hyp.g.moderator.volume.method.contrast[[3,7]],5))),"]")


primary.hyp.g.moderator.volume.method.plot=ggplot(data = primary.data.hyp%>%
                                                    mutate(set.rep.equated=ifelse(set.rep.equated=="set","Set Equated",
                                                                                  ifelse(set.rep.equated=="rep","Repetition Equated","Both Equated"))),
                                                  aes(x=avg.rir,
                                                      y=yi))+THEME+theme(legend.position = "right",
                                                                         plot.caption = element_text(face = "italic"))+
  facet_wrap(~set.rep.equated)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Method of Volume Equating x Estimated RIR",
       caption = primary.hyp.g.moderator.volume.method.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = primary.hyp.g.moderator.volume.method.emms%>%
                mutate(set.rep.equated=ifelse(set.rep.equated=="set","Set Equated",
                                              ifelse(set.rep.equated=="rep","Repetition Equated","Both Equated"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = primary.hyp.g.moderator.volume.method.emms%>%
                mutate(set.rep.equated=ifelse(set.rep.equated=="set","Set Equated",
                                              ifelse(set.rep.equated=="rep","Repetition Equated","Both Equated"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = primary.hyp.g.moderator.volume.method.emms%>%
              mutate(set.rep.equated=ifelse(set.rep.equated=="set","Set Equated",
                                            ifelse(set.rep.equated=="rep","Repetition Equated","Both Equated"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### g training status
primary.hyp.g.moderator.train.status=primary.hyp.g.best.model.prep%>%
  update(~. + avg.rir:train.status)%>%
  robust.rma.mv(cluster = primary.data.hyp$study)

primary.hyp.g.moderator.train.status.emms=primary.hyp.g.moderator.train.status%>%emmprep(at=list(avg.rir=0:23))%>%
  emmeans(~avg.rir+train.status,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.hyp.g.moderator.train.status)

primary.hyp.g.moderator.train.status.contrast=emmprep(primary.hyp.g.moderator.train.status,
                                                      at=list(avg.rir=primary.hyp.mean.rir))%>%
  emmeans(~avg.rir+train.status,
          weights = "prop")%>%
  contrast(interaction= c("consec","pairwise"))%>%
  confint()

primary.hyp.g.moderator.train.status.contrast.label=paste0('Interaction Contrast (Trained - Untrained) = ',
                                                           bquote(.(round(primary.hyp.g.moderator.train.status.contrast[[1,3]],5))),
                                                           " [95% CI: ",bquote(.(round(primary.hyp.g.moderator.train.status.contrast[[1,6]],5))),
                                                           ", ",bquote(.(round(primary.hyp.g.moderator.train.status.contrast[[1,7]],5))),"]")


primary.hyp.g.moderator.train.status.plot=ggplot(data = primary.data.hyp%>%
                                                   mutate(train.status=ifelse(train.status=="trained","Trained","Untrained")),
                                                 aes(x=avg.rir,
                                                     y=yi))+THEME+theme(legend.position = "right",
                                                                        plot.caption = element_text(face = "italic"))+
  facet_wrap(~train.status)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Training Status x Estimated RIR",
       caption = primary.hyp.g.moderator.train.status.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = primary.hyp.g.moderator.train.status.emms%>%
                mutate(train.status=ifelse(train.status=="trained","Trained","Untrained")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = primary.hyp.g.moderator.train.status.emms%>%
                mutate(train.status=ifelse(train.status=="trained","Trained","Untrained")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = primary.hyp.g.moderator.train.status.emms%>%
              mutate(train.status=ifelse(train.status=="trained","Trained","Untrained")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### g weeks
primary.hyp.g.moderator.weeks=primary.hyp.g.best.model.prep%>%
  update(~. + avg.rir:weeks)%>%
  robust.rma.mv(cluster = primary.data.hyp$study)

primary.hyp.g.moderator.weeks.emms=primary.hyp.g.moderator.weeks%>%
  emmprep(at=list(avg.rir=0:23,
                  weeks=c(4,8,12)))%>%
  emmeans(~avg.rir+weeks,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.hyp.g.moderator.weeks)

primary.hyp.g.moderator.weeks.contrast=emmprep(primary.hyp.g.moderator.weeks,
                                               at=list(avg.rir=primary.hyp.mean.rir,
                                                       weeks=primary.hyp.mean.weeks))%>%
  emmeans(~avg.rir+weeks,
          weights = "prop")%>%
  contrast("consec",interaction=TRUE)%>%
  confint()

primary.hyp.g.moderator.weeks.contrast.label=paste0('Interaction Contrast (Weeks) = ',
                                                    bquote(.(round(primary.hyp.g.moderator.weeks.contrast[[3]],5))),
                                                    " [95% CI: ",bquote(.(round(primary.hyp.g.moderator.weeks.contrast[[6]],5))),
                                                    ", ",bquote(.(round(primary.hyp.g.moderator.weeks.contrast[[7]],5))),"]")

primary.hyp.g.moderator.weeks.plot=ggplot(data = primary.data.hyp%>%
                                            mutate(weeks=ifelse(weeks<6,"4 Weeks",
                                                                ifelse(weeks>6&weeks<10,"8 Weeks","12 Weeks"))),
                                          aes(x=avg.rir,
                                              y=yi))+THEME+theme(plot.caption = element_text(face = "italic"))+
  facet_wrap(~fct_relevel(weeks,c("4 Weeks","8 Weeks","12 Weeks")))+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Intervention Duration (Weeks) x Estimated RIR",
       caption = primary.hyp.g.moderator.weeks.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = primary.hyp.g.moderator.weeks.emms%>%
                mutate(weeks=ifelse(weeks==4,"4 Weeks",
                                    ifelse(weeks==8,"8 Weeks","12 Weeks"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = primary.hyp.g.moderator.weeks.emms%>%
                mutate(weeks=ifelse(weeks==4,"4 Weeks",
                                    ifelse(weeks==8,"8 Weeks","12 Weeks"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = primary.hyp.g.moderator.weeks.emms%>%
              mutate(weeks=ifelse(weeks==4,"4 Weeks",
                                  ifelse(weeks==8,"8 Weeks","12 Weeks"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### g adjusted sets 
primary.hyp.g.moderator.adj.sets=primary.hyp.g.best.model.prep%>%
  update(~. + avg.rir*adj.sets.week)%>%
  robust.rma.mv(cluster = primary.data.hyp$study)

primary.hyp.g.moderator.adj.sets.emms=primary.hyp.g.moderator.adj.sets%>%
  emmprep(at=list(avg.rir=0:23,
                  adj.sets.week=c(8,16,24)))%>%
  emmeans(~avg.rir+adj.sets.week,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.hyp.g.moderator.adj.sets)

primary.hyp.g.moderator.adj.sets.contrast=emmprep(primary.hyp.g.moderator.adj.sets,
                                                  at=list(avg.rir=primary.hyp.mean.rir,
                                                          adj.sets.week=primary.hyp.mean.adj.sets.week))%>%
  emmeans(~avg.rir+adj.sets.week,
          weights = "prop")%>%
  contrast("consec",interaction=TRUE)%>%
  confint()

primary.hyp.g.moderator.adj.sets.contrast.label=paste0('Interaction Contrast (Adjusted Direct Sets) = ',
                                                       bquote(.(round(primary.hyp.g.moderator.adj.sets.contrast[[3]],5))),
                                                       " [95% CI: ",bquote(.(round(primary.hyp.g.moderator.adj.sets.contrast[[6]],5))),
                                                       ", ",bquote(.(round(primary.hyp.g.moderator.adj.sets.contrast[[7]],5))),"]")

primary.hyp.g.moderator.adj.sets.plot=ggplot(data = primary.data.hyp%>%
                                               mutate(adj.sets.week=ifelse(adj.sets.week<12,"8 Sets",
                                                                           ifelse(adj.sets.week>12&adj.sets.week<20,"16 Sets","24 Sets"))),
                                             aes(x=avg.rir,
                                                 y=yi))+THEME+theme(plot.caption = element_text(face = "italic"))+
  facet_wrap(~fct_relevel(adj.sets.week,c("8 Sets","16 Sets","24 Sets")))+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Adjusted Direct Sets Per Week Per Muscle Group x Estimated RIR",
       caption = primary.hyp.g.moderator.adj.sets.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = primary.hyp.g.moderator.adj.sets.emms%>%
                mutate(adj.sets.week=ifelse(adj.sets.week==8,"8 Sets",
                                            ifelse(adj.sets.week==16,"16 Sets","24 Sets"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = primary.hyp.g.moderator.adj.sets.emms%>%
                mutate(adj.sets.week=ifelse(adj.sets.week==8,"8 Sets",
                                            ifelse(adj.sets.week==16,"16 Sets","24 Sets"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = primary.hyp.g.moderator.adj.sets.emms%>%
              mutate(adj.sets.week=ifelse(adj.sets.week==8,"8 Sets",
                                          ifelse(adj.sets.week==16,"16 Sets","24 Sets"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### g frequency
primary.hyp.g.moderator.frequency.per.muscle=primary.hyp.g.best.model.prep%>%
  update(~. + avg.rir*frequency.per.muscle)%>%
  robust.rma.mv(cluster = primary.data.hyp$study)

primary.hyp.g.moderator.frequency.per.muscle.emms=primary.hyp.g.moderator.frequency.per.muscle%>%
  emmprep(at=list(avg.rir=0:23,
                  frequency.per.muscle=c(1,2,3)))%>%
  emmeans(~avg.rir+frequency.per.muscle,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.hyp.g.moderator.frequency.per.muscle)

primary.hyp.g.moderator.frequency.per.muscle.contrast=emmprep(primary.hyp.g.moderator.frequency.per.muscle,
                                                              at=list(avg.rir=primary.hyp.mean.rir,
                                                                      frequency.per.muscle=primary.hyp.mean.frequency.per.muscle))%>%
  emmeans(~avg.rir+frequency.per.muscle,
          weights = "prop")%>%
  contrast("consec",interaction=TRUE)%>%
  confint()

primary.hyp.g.moderator.frequency.per.muscle.contrast.label=paste0('Interaction Contrast (Direct Frequency) = ',
                                                                   bquote(.(round(primary.hyp.g.moderator.frequency.per.muscle.contrast[[3]],5))),
                                                                   " [95% CI: ",bquote(.(round(primary.hyp.g.moderator.frequency.per.muscle.contrast[[6]],5))),
                                                                   ", ",bquote(.(round(primary.hyp.g.moderator.frequency.per.muscle.contrast[[7]],5))),"]")

primary.hyp.g.moderator.frequency.per.muscle.plot=ggplot(data = primary.data.hyp%>%
                                                           mutate(frequency.per.muscle=ifelse(frequency.per.muscle<1.5,"1x Per Week",
                                                                                              ifelse(frequency.per.muscle>1.5&frequency.per.muscle<2.5,"2x Per Week","3x Per Week"))),
                                                         aes(x=avg.rir,
                                                             y=yi))+THEME+theme(plot.caption = element_text(face = "italic"))+
  facet_wrap(~frequency.per.muscle)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Direct Frequency Per Week Per Muscle Group x Estimated RIR",
       caption = primary.hyp.g.moderator.frequency.per.muscle.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = primary.hyp.g.moderator.frequency.per.muscle.emms%>%
                mutate(frequency.per.muscle=ifelse(frequency.per.muscle==1,"1x Per Week",
                                                   ifelse(frequency.per.muscle==2,"2x Per Week","3x Per Week"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = primary.hyp.g.moderator.frequency.per.muscle.emms%>%
                mutate(frequency.per.muscle=ifelse(frequency.per.muscle==1,"1x Per Week",
                                                   ifelse(frequency.per.muscle==2,"2x Per Week","3x Per Week"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = primary.hyp.g.moderator.frequency.per.muscle.emms%>%
              mutate(frequency.per.muscle=ifelse(frequency.per.muscle==1,"1x Per Week",
                                                 ifelse(frequency.per.muscle==2,"2x Per Week","3x Per Week"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### g age
primary.hyp.g.moderator.age=primary.hyp.g.best.model.prep%>%
  update(~. + avg.rir*age)%>%
  robust.rma.mv(cluster = primary.data.hyp$study)

primary.hyp.g.moderator.age.emms=primary.hyp.g.moderator.age%>%
  emmprep(at=list(avg.rir=0:23,
                  age=c(20,40,60)))%>%
  emmeans(~avg.rir+age,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.hyp.g.moderator.age)

primary.hyp.g.moderator.age.contrast=emmprep(primary.hyp.g.moderator.age,
                                             at=list(avg.rir=primary.hyp.mean.rir,
                                                     age=primary.hyp.mean.age))%>%
  emmeans(~avg.rir+age,
          weights = "prop")%>%
  contrast("consec",interaction=TRUE)%>%
  confint()

primary.hyp.g.moderator.age.contrast.label=paste0('Interaction Contrast (Age) = ',
                                                  bquote(.(round(primary.hyp.g.moderator.age.contrast[[3]],5))),
                                                  " [95% CI: ",bquote(.(round(primary.hyp.g.moderator.age.contrast[[6]],5))),
                                                  ", ",bquote(.(round(primary.hyp.g.moderator.age.contrast[[7]],5))),"]")

primary.hyp.g.moderator.age.plot=ggplot(data = primary.data.hyp%>%
                                          mutate(age=ifelse(age<30,"20 Years Old",
                                                            ifelse(age>30&age<50,"40 Years Old","60 Years Old"))),
                                        aes(x=avg.rir,
                                            y=yi))+THEME+theme(plot.caption = element_text(face = "italic"))+
  facet_wrap(~age)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Age x Estimated RIR",
       caption = primary.hyp.g.moderator.age.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = primary.hyp.g.moderator.age.emms%>%
                mutate(age=ifelse(age==20,"20 Years Old",
                                  ifelse(age==40,"40 Years Old","60 Years Old"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = primary.hyp.g.moderator.age.emms%>%
                mutate(age=ifelse(age==20,"20 Years Old",
                                  ifelse(age==40,"40 Years Old","60 Years Old"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = primary.hyp.g.moderator.age.emms%>%
              mutate(age=ifelse(age==20,"20 Years Old",
                                ifelse(age==40,"40 Years Old","60 Years Old"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### g sex
primary.hyp.g.moderator.sex.percent.male=primary.hyp.g.best.model.prep%>%
  update(~. + avg.rir*sex.percent.male)%>%
  robust.rma.mv(cluster = primary.data.hyp$study)

primary.hyp.g.moderator.sex.percent.male.emms=primary.hyp.g.moderator.sex.percent.male%>%
  emmprep(at=list(avg.rir=0:23,
                  sex.percent.male=c(0,100)))%>%
  emmeans(~avg.rir+sex.percent.male,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.hyp.g.moderator.sex.percent.male)

primary.hyp.g.moderator.sex.percent.male.contrast=emmprep(primary.hyp.g.moderator.sex.percent.male,
                                                          at=list(avg.rir=primary.hyp.mean.rir,
                                                                  sex.percent.male=primary.hyp.mean.sex))%>%
  emmeans(~avg.rir+sex.percent.male,
          weights = "prop")%>%
  contrast("consec",interaction=TRUE)%>%
  confint()

primary.hyp.g.moderator.sex.percent.male.contrast.label=paste0('Interaction Contrast (Sex) = ',
                                                               bquote(.(round(primary.hyp.g.moderator.sex.percent.male.contrast[[3]],5))),
                                                               " [95% CI: ",bquote(.(round(primary.hyp.g.moderator.sex.percent.male.contrast[[6]],5))),
                                                               ", ",bquote(.(round(primary.hyp.g.moderator.sex.percent.male.contrast[[7]],5))),"]")

primary.hyp.g.moderator.sex.percent.male.plot=ggplot(data = primary.data.hyp%>%
                                                       mutate(sex.percent.male=ifelse(sex.percent.male<50,"Female","Male")),
                                                     aes(x=avg.rir,
                                                         y=yi))+THEME+theme(plot.caption = element_text(face = "italic"))+
  facet_wrap(~sex.percent.male)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Sex (% Male) x Estimated RIR",
       caption = primary.hyp.g.moderator.sex.percent.male.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = primary.hyp.g.moderator.sex.percent.male.emms%>%
                mutate(sex.percent.male=ifelse(sex.percent.male==0,"Female","Male")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = primary.hyp.g.moderator.sex.percent.male.emms%>%
                mutate(sex.percent.male=ifelse(sex.percent.male==0,"Female","Male")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = primary.hyp.g.moderator.sex.percent.male.emms%>%
              mutate(sex.percent.male=ifelse(sex.percent.male==0,"Female","Male")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### g study design
primary.hyp.g.moderator.within.between.design=primary.hyp.g.best.model.prep%>%
  update(~. + avg.rir*within.between.design)%>%
  robust.rma.mv(cluster = primary.data.hyp$study)

primary.hyp.g.moderator.within.between.design.emms=primary.hyp.g.moderator.within.between.design%>%emmprep(at=list(avg.rir=0:23))%>%
  emmeans(~avg.rir+within.between.design,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.hyp.g.moderator.within.between.design)

primary.hyp.g.moderator.within.between.design.contrast=emmprep(primary.hyp.g.moderator.within.between.design,
                                                               at=list(avg.rir=primary.hyp.mean.rir))%>%
  emmeans(~avg.rir+within.between.design,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()

primary.hyp.g.moderator.within.between.design.contrast.label=paste0('Interaction Contrast (Between - Within) = ',
                                                                    bquote(.(round(primary.hyp.g.moderator.within.between.design.contrast[[1,3]],5))),
                                                                    " [95% CI: ",bquote(.(round(primary.hyp.g.moderator.within.between.design.contrast[[1,6]],5))),
                                                                    ", ",bquote(.(round(primary.hyp.g.moderator.within.between.design.contrast[[1,7]],5))),"]")


primary.hyp.g.moderator.within.between.design.plot=ggplot(data = primary.data.hyp%>%
                                                            mutate(within.between.design=ifelse(within.between.design=="Between","Between Subjects","Within Subjects")),
                                                          aes(x=avg.rir,
                                                              y=yi))+THEME+theme(legend.position = "right",
                                                                                 plot.caption = element_text(face = "italic"))+
  facet_wrap(~within.between.design)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Study Design x Estimated RIR",
       caption = primary.hyp.g.moderator.within.between.design.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = primary.hyp.g.moderator.within.between.design.emms%>%
                mutate(within.between.design=ifelse(within.between.design=="Between","Between Subjects","Within Subjects")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = primary.hyp.g.moderator.within.between.design.emms%>%
                mutate(within.between.design=ifelse(within.between.design=="Between","Between Subjects","Within Subjects")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = primary.hyp.g.moderator.within.between.design.emms%>%
              mutate(within.between.design=ifelse(within.between.design=="Between","Between Subjects","Within Subjects")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### g upper / lower outcome
primary.hyp.g.moderator.upper.lower.other=primary.hyp.g.best.model.prep%>%
  update(~. + avg.rir*upper.lower.other)%>%
  robust.rma.mv(cluster = primary.data.hyp$study)

primary.hyp.g.moderator.upper.lower.other.emms=primary.hyp.g.moderator.upper.lower.other%>%emmprep(at=list(avg.rir=0:23))%>%
  emmeans(~avg.rir+upper.lower.other,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.hyp.g.moderator.upper.lower.other)

primary.hyp.g.moderator.upper.lower.other.contrast=emmprep(primary.hyp.g.moderator.upper.lower.other,
                                                           at=list(avg.rir=primary.hyp.mean.rir))%>%
  emmeans(~avg.rir+upper.lower.other,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()

primary.hyp.g.moderator.upper.lower.other.contrast.label=paste0('Interaction Contrast (Lower - Upper) = ',
                                                                bquote(.(round(primary.hyp.g.moderator.upper.lower.other.contrast[[1,3]],5))),
                                                                " [95% CI: ",bquote(.(round(primary.hyp.g.moderator.upper.lower.other.contrast[[1,6]],5))),
                                                                ", ",bquote(.(round(primary.hyp.g.moderator.upper.lower.other.contrast[[1,7]],5))),"]")


primary.hyp.g.moderator.upper.lower.other.plot=ggplot(data = primary.data.hyp%>%
                                                        mutate(upper.lower.other=ifelse(upper.lower.other=="Lower","Lower Body Outcomes","Upper Body Outcomes")),
                                                      aes(x=avg.rir,
                                                          y=yi))+THEME+theme(legend.position = "right",
                                                                             plot.caption = element_text(face = "italic"))+
  facet_wrap(~upper.lower.other)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Upper/Lower Body Outcome x Estimated RIR",
       caption = primary.hyp.g.moderator.upper.lower.other.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = primary.hyp.g.moderator.upper.lower.other.emms%>%
                mutate(upper.lower.other=ifelse(upper.lower.other=="Lower","Lower Body Outcomes","Upper Body Outcomes")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = primary.hyp.g.moderator.upper.lower.other.emms%>%
                mutate(upper.lower.other=ifelse(upper.lower.other=="Lower","Lower Body Outcomes","Upper Body Outcomes")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = primary.hyp.g.moderator.upper.lower.other.emms%>%
              mutate(upper.lower.other=ifelse(upper.lower.other=="Lower","Lower Body Outcomes","Upper Body Outcomes")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### g exercise selection
primary.hyp.g.moderator.train.exercise=primary.hyp.g.best.model.prep%>%
  update(~. + avg.rir*train.exercise)%>%
  robust.rma.mv(cluster = primary.data.hyp$study)

primary.hyp.g.moderator.train.exercise.emms=primary.hyp.g.moderator.train.exercise%>%emmprep(at=list(avg.rir=0:23))%>%
  emmeans(~avg.rir+train.exercise,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.hyp.g.moderator.train.exercise)

primary.hyp.g.moderator.train.exercise.contrast=emmprep(primary.hyp.g.moderator.train.exercise,
                                                        at=list(avg.rir=primary.hyp.mean.rir))%>%
  emmeans(~avg.rir+train.exercise,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()

primary.hyp.g.moderator.train.exercise.contrast.label=paste0('Interaction Contrast (Both - Multi-joint) = ',
                                                             bquote(.(round(primary.hyp.g.moderator.train.exercise.contrast[[1,3]],5))),
                                                             " [95% CI: ",bquote(.(round(primary.hyp.g.moderator.train.exercise.contrast[[1,6]],5))),
                                                             ", ",bquote(.(round(primary.hyp.g.moderator.train.exercise.contrast[[1,7]],5))),"]\n",
                                                             'Interaction Contrast (Both - Single-joint) = ',
                                                             bquote(.(round(primary.hyp.g.moderator.train.exercise.contrast[[2,3]],5))),
                                                             " [95% CI: ",bquote(.(round(primary.hyp.g.moderator.train.exercise.contrast[[2,6]],5))),
                                                             ", ",bquote(.(round(primary.hyp.g.moderator.train.exercise.contrast[[2,7]],5))),"]\n",
                                                             'Interaction Contrast (Multi-joint - Single-joint) = ',
                                                             bquote(.(round(primary.hyp.g.moderator.train.exercise.contrast[[3,3]],5))),
                                                             " [95% CI: ",bquote(.(round(primary.hyp.g.moderator.train.exercise.contrast[[3,6]],5))),
                                                             ", ",bquote(.(round(primary.hyp.g.moderator.train.exercise.contrast[[3,7]],5))),"]")


primary.hyp.g.moderator.train.exercise.plot=ggplot(data = primary.data.hyp%>%
                                                     mutate(train.exercise=ifelse(train.exercise=="Both","Both",
                                                                                  ifelse(train.exercise=="Multi-joint","Multi-joint Exercises","Single-joint Exercises"))),
                                                   aes(x=avg.rir,
                                                       y=yi))+THEME+theme(legend.position = "right",
                                                                          plot.caption = element_text(face = "italic"))+
  facet_wrap(~train.exercise)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Exercise Selection x Estimated RIR",
       caption = primary.hyp.g.moderator.train.exercise.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = primary.hyp.g.moderator.train.exercise.emms%>%
                mutate(train.exercise=ifelse(train.exercise=="Both","Both",
                                             ifelse(train.exercise=="Multi-joint","Multi-joint Exercises","Single-joint Exercises"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = primary.hyp.g.moderator.train.exercise.emms%>%
                mutate(train.exercise=ifelse(train.exercise=="Both","Both",
                                             ifelse(train.exercise=="Multi-joint","Multi-joint Exercises","Single-joint Exercises"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = primary.hyp.g.moderator.train.exercise.emms%>%
              mutate(train.exercise=ifelse(train.exercise=="Both","Both",
                                           ifelse(train.exercise=="Multi-joint","Multi-joint Exercises","Single-joint Exercises"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### g concurrent training
primary.hyp.g.moderator.formal.cardio.intervention=primary.hyp.g.best.model.prep%>%
  update(~. + avg.rir*formal.cardio.intervention)%>%
  robust.rma.mv(cluster = primary.data.hyp$study)

primary.hyp.g.moderator.formal.cardio.intervention.emms=primary.hyp.g.moderator.formal.cardio.intervention%>%emmprep(at=list(avg.rir=0:23))%>%
  emmeans(~avg.rir+formal.cardio.intervention,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.hyp.g.moderator.formal.cardio.intervention)

primary.hyp.g.moderator.formal.cardio.intervention.contrast=emmprep(primary.hyp.g.moderator.formal.cardio.intervention,
                                                                    at=list(avg.rir=primary.hyp.mean.rir))%>%
  emmeans(~avg.rir+formal.cardio.intervention,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()

primary.hyp.g.moderator.formal.cardio.intervention.contrast.label=paste0('Interaction Contrast (RT Only - Concurrent) = ',
                                                                         bquote(.(round(primary.hyp.g.moderator.formal.cardio.intervention.contrast[[1,3]],5))),
                                                                         " [95% CI: ",bquote(.(round(primary.hyp.g.moderator.formal.cardio.intervention.contrast[[1,6]],5))),
                                                                         ", ",bquote(.(round(primary.hyp.g.moderator.formal.cardio.intervention.contrast[[1,7]],5))),"]")


primary.hyp.g.moderator.formal.cardio.intervention.plot=ggplot(data = primary.data.hyp%>%
                                                                 mutate(formal.cardio.intervention=ifelse(formal.cardio.intervention=="No","Only Resistance Training","Concurrent Training")),
                                                               aes(x=avg.rir,
                                                                   y=yi))+THEME+theme(legend.position = "right",
                                                                                      plot.caption = element_text(face = "italic"))+
  facet_wrap(~formal.cardio.intervention)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Concurrent Training x Estimated RIR",
       caption = primary.hyp.g.moderator.formal.cardio.intervention.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = primary.hyp.g.moderator.formal.cardio.intervention.emms%>%
                mutate(formal.cardio.intervention=ifelse(formal.cardio.intervention=="No","Only Resistance Training","Concurrent Training")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = primary.hyp.g.moderator.formal.cardio.intervention.emms%>%
                mutate(formal.cardio.intervention=ifelse(formal.cardio.intervention=="No","Only Resistance Training","Concurrent Training")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = primary.hyp.g.moderator.formal.cardio.intervention.emms%>%
              mutate(formal.cardio.intervention=ifelse(formal.cardio.intervention=="No","Only Resistance Training","Concurrent Training")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### g progressive overload
primary.hyp.g.moderator.progression=primary.hyp.g.best.model.prep%>%
  update(~. + avg.rir*progression)%>%
  robust.rma.mv(cluster = primary.data.hyp$study)

primary.hyp.g.moderator.progression.emms=primary.hyp.g.moderator.progression%>%emmprep(at=list(avg.rir=0:23))%>%
  emmeans(~avg.rir+progression,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.hyp.g.moderator.progression)

primary.hyp.g.moderator.progression.contrast=emmprep(primary.hyp.g.moderator.progression,
                                                     at=list(avg.rir=primary.hyp.mean.rir))%>%
  emmeans(~avg.rir+progression,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()

primary.hyp.g.moderator.progression.contrast.label=paste0('Interaction Contrast (No Overload - Progressive Overload) = ',
                                                          bquote(.(round(primary.hyp.g.moderator.progression.contrast[[1,3]],5))),
                                                          " [95% CI: ",bquote(.(round(primary.hyp.g.moderator.progression.contrast[[1,6]],5))),
                                                          ", ",bquote(.(round(primary.hyp.g.moderator.progression.contrast[[1,7]],5))),"]")


primary.hyp.g.moderator.progression.plot=ggplot(data = primary.data.hyp%>%
                                                  mutate(progression=ifelse(progression=="No","No Progressive Overload","Progressive Overload")),
                                                aes(x=avg.rir,
                                                    y=yi))+THEME+theme(legend.position = "right",
                                                                       plot.caption = element_text(face = "italic"))+
  facet_wrap(~progression)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Progressive Overload x Estimated RIR",
       caption = primary.hyp.g.moderator.progression.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = primary.hyp.g.moderator.progression.emms%>%
                mutate(progression=ifelse(progression=="No","No Progressive Overload","Progressive Overload")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = primary.hyp.g.moderator.progression.emms%>%
                mutate(progression=ifelse(progression=="No","No Progressive Overload","Progressive Overload")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = primary.hyp.g.moderator.progression.emms%>%
              mutate(progression=ifelse(progression=="No","No Progressive Overload","Progressive Overload")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### g failure definition
primary.hyp.g.moderator.failure.definition=primary.hyp.g.best.model.prep%>%
  update(~. + avg.rir*failure.definition)%>%
  robust.rma.mv(cluster = primary.data.hyp$study)

primary.hyp.g.moderator.failure.definition.emms=primary.hyp.g.moderator.failure.definition%>%emmprep(at=list(avg.rir=0:23))%>%
  emmeans(~avg.rir+failure.definition,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.hyp.g.moderator.failure.definition)

primary.hyp.g.moderator.failure.definition.contrast=emmprep(primary.hyp.g.moderator.failure.definition,
                                                            at=list(avg.rir=primary.hyp.mean.rir))%>%
  emmeans(~avg.rir+failure.definition,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()

primary.hyp.g.moderator.failure.definition.contrast.label=paste0('Interaction Contrast (No Definition - Definition Provided) = ',
                                                                 bquote(.(round(primary.hyp.g.moderator.failure.definition.contrast[[1,3]],5))),
                                                                 " [95% CI: ",bquote(.(round(primary.hyp.g.moderator.failure.definition.contrast[[1,6]],5))),
                                                                 ", ",bquote(.(round(primary.hyp.g.moderator.failure.definition.contrast[[1,7]],5))),"]")


primary.hyp.g.moderator.failure.definition.plot=ggplot(data = primary.data.hyp%>%
                                                         mutate(failure.definition=ifelse(failure.definition=="No","No Failure Definition","Failure Definition Provided")),
                                                       aes(x=avg.rir,
                                                           y=yi))+THEME+theme(legend.position = "right",
                                                                              plot.caption = element_text(face = "italic"))+
  facet_wrap(~failure.definition)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Failure Definition x Estimated RIR",
       caption = primary.hyp.g.moderator.failure.definition.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = primary.hyp.g.moderator.failure.definition.emms%>%
                mutate(failure.definition=ifelse(failure.definition=="No","No Failure Definition","Failure Definition Provided")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = primary.hyp.g.moderator.failure.definition.emms%>%
                mutate(failure.definition=ifelse(failure.definition=="No","No Failure Definition","Failure Definition Provided")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = primary.hyp.g.moderator.failure.definition.emms%>%
              mutate(failure.definition=ifelse(failure.definition=="No","No Failure Definition","Failure Definition Provided")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### g set structure
primary.hyp.g.moderator.alternative.set.structure=primary.hyp.g.best.model.prep%>%
  update(~. + avg.rir*alternative.set.structure)%>%
  robust.rma.mv(cluster = primary.data.hyp$study)

primary.hyp.g.moderator.alternative.set.structure.emms=primary.hyp.g.moderator.alternative.set.structure%>%emmprep(at=list(avg.rir=0:23))%>%
  emmeans(~avg.rir+alternative.set.structure,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.hyp.g.moderator.alternative.set.structure)

primary.hyp.g.moderator.alternative.set.structure.contrast=emmprep(primary.hyp.g.moderator.alternative.set.structure,
                                                                   at=list(avg.rir=primary.hyp.mean.rir))%>%
  emmeans(~avg.rir+alternative.set.structure,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()

primary.hyp.g.moderator.alternative.set.structure.contrast.label=paste0('Interaction Contrast (Traditional - Alternative) = ',
                                                                        bquote(.(round(primary.hyp.g.moderator.alternative.set.structure.contrast[[1,3]],5))),
                                                                        " [95% CI: ",bquote(.(round(primary.hyp.g.moderator.alternative.set.structure.contrast[[1,6]],5))),
                                                                        ", ",bquote(.(round(primary.hyp.g.moderator.alternative.set.structure.contrast[[1,7]],5))),"]")


primary.hyp.g.moderator.alternative.set.structure.plot=ggplot(data = primary.data.hyp%>%
                                                                mutate(alternative.set.structure=ifelse(alternative.set.structure=="No","Traditional Set Structure","Alternative Set Structure")),
                                                              aes(x=avg.rir,
                                                                  y=yi))+THEME+theme(legend.position = "right",
                                                                                     plot.caption = element_text(face = "italic"))+
  facet_wrap(~alternative.set.structure)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Set Structure x Estimated RIR",
       caption = primary.hyp.g.moderator.alternative.set.structure.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = primary.hyp.g.moderator.alternative.set.structure.emms%>%
                mutate(alternative.set.structure=ifelse(alternative.set.structure=="No","Traditional Set Structure","Alternative Set Structure")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = primary.hyp.g.moderator.alternative.set.structure.emms%>%
                mutate(alternative.set.structure=ifelse(alternative.set.structure=="No","Traditional Set Structure","Alternative Set Structure")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = primary.hyp.g.moderator.alternative.set.structure.emms%>%
              mutate(alternative.set.structure=ifelse(alternative.set.structure=="No","Traditional Set Structure","Alternative Set Structure")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### g intended velocity
primary.hyp.g.moderator.train.intent=primary.hyp.g.best.model.prep%>%
  update(~. + avg.rir*train.intent)%>%
  robust.rma.mv(cluster = primary.data.hyp$study)

primary.hyp.g.moderator.train.intent.emms=primary.hyp.g.moderator.train.intent%>%emmprep(at=list(avg.rir=0:23))%>%
  emmeans(~avg.rir+train.intent,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.hyp.g.moderator.train.intent)

primary.hyp.g.moderator.train.intent.contrast=emmprep(primary.hyp.g.moderator.train.intent,
                                                      at=list(avg.rir=primary.hyp.mean.rir))%>%
  emmeans(~avg.rir+train.intent,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()

primary.hyp.g.moderator.train.intent.contrast.label=paste0('Interaction Contrast (Submaximal - Maximal) = ',
                                                           bquote(.(round(primary.hyp.g.moderator.train.intent.contrast[[1,3]],5))),
                                                           " [95% CI: ",bquote(.(round(primary.hyp.g.moderator.train.intent.contrast[[1,6]],5))),
                                                           ", ",bquote(.(round(primary.hyp.g.moderator.train.intent.contrast[[1,7]],5))),"]")


primary.hyp.g.moderator.train.intent.plot=ggplot(data = primary.data.hyp%>%
                                                   mutate(train.intent=ifelse(train.intent=="No","Submaximal Intended Concentric Velocity","Maximal Intended Concentric Velocity")),
                                                 aes(x=avg.rir,
                                                     y=yi))+THEME+theme(legend.position = "right",
                                                                        plot.caption = element_text(face = "italic"))+
  facet_wrap(~train.intent)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Intended Concentric Velocity x Estimated RIR",
       caption = primary.hyp.g.moderator.train.intent.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = primary.hyp.g.moderator.train.intent.emms%>%
                mutate(train.intent=ifelse(train.intent=="No","Submaximal Intended Concentric Velocity","Maximal Intended Concentric Velocity")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = primary.hyp.g.moderator.train.intent.emms%>%
                mutate(train.intent=ifelse(train.intent=="No","Submaximal Intended Concentric Velocity","Maximal Intended Concentric Velocity")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = primary.hyp.g.moderator.train.intent.emms%>%
              mutate(train.intent=ifelse(train.intent=="No","Submaximal Intended Concentric Velocity","Maximal Intended Concentric Velocity")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### combine g plots
primary.hyp.g.moderator.doc=primary.hyp.g.moderator.load.plot/plot_spacer()/primary.hyp.g.moderator.volume.method.plot/plot_spacer()/
  primary.hyp.g.moderator.train.status.plot/plot_spacer()/primary.hyp.g.moderator.weeks.plot/plot_spacer()/
  primary.hyp.g.moderator.adj.sets.plot/plot_spacer()/primary.hyp.g.moderator.frequency.per.muscle.plot/plot_spacer()/
  primary.hyp.g.moderator.age.plot/plot_spacer()/primary.hyp.g.moderator.sex.percent.male.plot/plot_spacer()/
  primary.hyp.g.moderator.within.between.design.plot/plot_spacer()/primary.hyp.g.moderator.upper.lower.other.plot/plot_spacer()/
  primary.hyp.g.moderator.train.exercise.plot/plot_spacer()/primary.hyp.g.moderator.formal.cardio.intervention.plot/plot_spacer()/
  primary.hyp.g.moderator.progression.plot/plot_spacer()/primary.hyp.g.moderator.failure.definition.plot/plot_spacer()/
  primary.hyp.g.moderator.alternative.set.structure.plot/plot_spacer()/primary.hyp.g.moderator.train.intent.plot&
  plot_annotation(title = "Moderator Analyses from Linear Multilevel Model for Muscle Hypertrophy (Standardized Mean Change)")&
  THEME

ggsave(primary.hyp.g.moderator.doc,device = "pdf",filename="primary.hyp.g.moderator.doc.pdf",width = 12,height = 136.5,limitsize = FALSE,
       path = "/Users/zacrobinson/Documents/For Publication PTF Meta/Plots")


#### rom load
primary.hyp.rom.moderator.load=primary.hyp.rom.best.model.prep%>%
  update(~. + avg.rir:load.set)%>%
  robust.rma.mv(cluster = primary.data.hyp$study)

primary.hyp.rom.moderator.load.emms=primary.hyp.rom.moderator.load%>%
  emmprep(at=list(avg.rir=0:23,
                  load.set=c(30,60,90)))%>%
  emmeans(~avg.rir+load.set,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.hyp.rom.moderator.load)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

primary.hyp.rom.moderator.load.contrast=emmprep(primary.hyp.rom.moderator.load,
                                                at=list(avg.rir=primary.hyp.mean.rir,
                                                        load.set=primary.hyp.mean.load.set))%>%
  emmeans(~avg.rir+load.set,
          weights = "prop")%>%
  contrast("consec",interaction=TRUE)%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

primary.hyp.rom.moderator.load.contrast.label=paste0('Interaction Contrast (Load) = ',
                                                     bquote(.(round(primary.hyp.rom.moderator.load.contrast[[3]],5))),
                                                     " [95% CI: ",bquote(.(round(primary.hyp.rom.moderator.load.contrast[[6]],5))),
                                                     ", ",bquote(.(round(primary.hyp.rom.moderator.load.contrast[[7]],5))),"]")

primary.hyp.rom.moderator.load.plot=ggplot(data = primary.data.hyp%>%
                                             mutate(load.set=ifelse(load.set<45,"30% of 1RM",
                                                                    ifelse(load.set>45&load.set<75,"60% of 1RM","90% of 1RM"))),
                                           aes(x=avg.rir,
                                               y=yi.rom.exp))+THEME+theme(plot.caption = element_text(face = "italic"))+
  facet_wrap(~load.set)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Load (% of 1RM) x Estimated RIR",
       caption = primary.hyp.rom.moderator.load.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = primary.hyp.rom.moderator.load.emms%>%
                mutate(load.set=ifelse(load.set==30,"30% of 1RM",
                                       ifelse(load.set==60,"60% of 1RM","90% of 1RM"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = primary.hyp.rom.moderator.load.emms%>%
                mutate(load.set=ifelse(load.set==30,"30% of 1RM",
                                       ifelse(load.set==60,"60% of 1RM","90% of 1RM"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = primary.hyp.rom.moderator.load.emms%>%
              mutate(load.set=ifelse(load.set==30,"30% of 1RM",
                                     ifelse(load.set==60,"60% of 1RM","90% of 1RM"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### rom volume equating method
primary.hyp.rom.moderator.volume.method=primary.hyp.rom.best.model.prep%>%
  update(~. + avg.rir:set.rep.equated)%>%
  robust.rma.mv(cluster = primary.data.hyp$study)

primary.hyp.rom.moderator.volume.method.emms=primary.hyp.rom.moderator.volume.method%>%emmprep(at=list(avg.rir=0:23))%>%
  emmeans(~avg.rir+set.rep.equated,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.hyp.rom.moderator.volume.method)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

primary.hyp.rom.moderator.volume.method.contrast=emmprep(primary.hyp.rom.moderator.volume.method,
                                                         at=list(avg.rir=primary.hyp.mean.rir))%>%
  emmeans(~avg.rir+set.rep.equated,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

primary.hyp.rom.moderator.volume.method.contrast.label=paste0('Interaction Contrast (Set - Rep) = ',
                                                              bquote(.(round(primary.hyp.rom.moderator.volume.method.contrast[[1,3]],5))),
                                                              " [95% CI: ",bquote(.(round(primary.hyp.rom.moderator.volume.method.contrast[[1,6]],5))),
                                                              ", ",bquote(.(round(primary.hyp.rom.moderator.volume.method.contrast[[1,7]],5))),"]\n",
                                                              'Interaction Contrast (Set - Both) = ',
                                                              bquote(.(round(primary.hyp.rom.moderator.volume.method.contrast[[2,3]],5))),
                                                              " [95% CI: ",bquote(.(round(primary.hyp.rom.moderator.volume.method.contrast[[2,6]],5))),
                                                              ", ",bquote(.(round(primary.hyp.rom.moderator.volume.method.contrast[[2,7]],5))),"]\n",
                                                              'Interaction Contrast (Rep - Both) = ',
                                                              bquote(.(round(primary.hyp.rom.moderator.volume.method.contrast[[3,3]],5))),
                                                              " [95% CI: ",bquote(.(round(primary.hyp.rom.moderator.volume.method.contrast[[3,6]],5))),
                                                              ", ",bquote(.(round(primary.hyp.rom.moderator.volume.method.contrast[[3,7]],5))),"]")


primary.hyp.rom.moderator.volume.method.plot=ggplot(data = primary.data.hyp%>%
                                                      mutate(set.rep.equated=ifelse(set.rep.equated=="set","Set Equated",
                                                                                    ifelse(set.rep.equated=="rep","Repetition Equated","Both Equated"))),
                                                    aes(x=avg.rir,
                                                        y=yi.rom.exp))+THEME+theme(legend.position = "right",
                                                                                   plot.caption = element_text(face = "italic"))+
  facet_wrap(~set.rep.equated)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Method of Volume Equating x Estimated RIR",
       caption = primary.hyp.rom.moderator.volume.method.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = primary.hyp.rom.moderator.volume.method.emms%>%
                mutate(set.rep.equated=ifelse(set.rep.equated=="set","Set Equated",
                                              ifelse(set.rep.equated=="rep","Repetition Equated","Both Equated"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = primary.hyp.rom.moderator.volume.method.emms%>%
                mutate(set.rep.equated=ifelse(set.rep.equated=="set","Set Equated",
                                              ifelse(set.rep.equated=="rep","Repetition Equated","Both Equated"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = primary.hyp.rom.moderator.volume.method.emms%>%
              mutate(set.rep.equated=ifelse(set.rep.equated=="set","Set Equated",
                                            ifelse(set.rep.equated=="rep","Repetition Equated","Both Equated"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### rom training status
primary.hyp.rom.moderator.train.status=primary.hyp.rom.best.model.prep%>%
  update(~. + avg.rir:train.status)%>%
  robust.rma.mv(cluster = primary.data.hyp$study)

primary.hyp.rom.moderator.train.status.emms=primary.hyp.rom.moderator.train.status%>%emmprep(at=list(avg.rir=0:23))%>%
  emmeans(~avg.rir+train.status,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.hyp.rom.moderator.train.status)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

primary.hyp.rom.moderator.train.status.contrast=emmprep(primary.hyp.rom.moderator.train.status,
                                                        at=list(avg.rir=primary.hyp.mean.rir))%>%
  emmeans(~avg.rir+train.status,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

primary.hyp.rom.moderator.train.status.contrast.label=paste0('Interaction Contrast (Trained - Untrained) = ',
                                                             bquote(.(round(primary.hyp.rom.moderator.train.status.contrast[[1,3]],5))),
                                                             " [95% CI: ",bquote(.(round(primary.hyp.rom.moderator.train.status.contrast[[1,6]],5))),
                                                             ", ",bquote(.(round(primary.hyp.rom.moderator.train.status.contrast[[1,7]],5))),"]")


primary.hyp.rom.moderator.train.status.plot=ggplot(data = primary.data.hyp%>%
                                                     mutate(train.status=ifelse(train.status=="trained","Trained","Untrained")),
                                                   aes(x=avg.rir,
                                                       y=yi.rom.exp))+THEME+theme(legend.position = "right",
                                                                                  plot.caption = element_text(face = "italic"))+
  facet_wrap(~train.status)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Training Status x Estimated RIR",
       caption = primary.hyp.rom.moderator.train.status.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = primary.hyp.rom.moderator.train.status.emms%>%
                mutate(train.status=ifelse(train.status=="trained","Trained","Untrained")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = primary.hyp.rom.moderator.train.status.emms%>%
                mutate(train.status=ifelse(train.status=="trained","Trained","Untrained")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = primary.hyp.rom.moderator.train.status.emms%>%
              mutate(train.status=ifelse(train.status=="trained","Trained","Untrained")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### rom weeks
primary.hyp.rom.moderator.weeks=primary.hyp.rom.best.model.prep%>%
  update(~. + avg.rir:weeks)%>%
  robust.rma.mv(cluster = primary.data.hyp$study)

primary.hyp.rom.moderator.weeks.emms=primary.hyp.rom.moderator.weeks%>%
  emmprep(at=list(avg.rir=0:23,
                  weeks=c(4,8,12)))%>%
  emmeans(~avg.rir+weeks,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.hyp.rom.moderator.weeks)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

primary.hyp.rom.moderator.weeks.contrast=emmprep(primary.hyp.rom.moderator.weeks,
                                                 at=list(avg.rir=primary.hyp.mean.rir,
                                                         weeks=primary.hyp.mean.weeks))%>%
  emmeans(~avg.rir+weeks,
          weights = "prop")%>%
  contrast("consec",interaction=TRUE)%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

primary.hyp.rom.moderator.weeks.contrast.label=paste0('Interaction Contrast (Weeks) = ',
                                                      bquote(.(round(primary.hyp.rom.moderator.weeks.contrast[[3]],5))),
                                                      " [95% CI: ",bquote(.(round(primary.hyp.rom.moderator.weeks.contrast[[6]],5))),
                                                      ", ",bquote(.(round(primary.hyp.rom.moderator.weeks.contrast[[7]],5))),"]")

primary.hyp.rom.moderator.weeks.plot=ggplot(data = primary.data.hyp%>%
                                              mutate(weeks=ifelse(weeks<6,"4 Weeks",
                                                                  ifelse(weeks>6&weeks<10,"8 Weeks","12 Weeks"))),
                                            aes(x=avg.rir,
                                                y=yi.rom.exp))+THEME+theme(plot.caption = element_text(face = "italic"))+
  facet_wrap(~fct_relevel(weeks,c("4 Weeks","8 Weeks","12 Weeks")))+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Intervention Duration (Weeks) x Estimated RIR",
       caption = primary.hyp.rom.moderator.weeks.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = primary.hyp.rom.moderator.weeks.emms%>%
                mutate(weeks=ifelse(weeks==4,"4 Weeks",
                                    ifelse(weeks==8,"8 Weeks","12 Weeks"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = primary.hyp.rom.moderator.weeks.emms%>%
                mutate(weeks=ifelse(weeks==4,"4 Weeks",
                                    ifelse(weeks==8,"8 Weeks","12 Weeks"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = primary.hyp.rom.moderator.weeks.emms%>%
              mutate(weeks=ifelse(weeks==4,"4 Weeks",
                                  ifelse(weeks==8,"8 Weeks","12 Weeks"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### rom adjusted sets 
primary.hyp.rom.moderator.adj.sets=primary.hyp.rom.best.model.prep%>%
  update(~. + avg.rir*adj.sets.week)%>%
  robust.rma.mv(cluster = primary.data.hyp$study)

primary.hyp.rom.moderator.adj.sets.emms=primary.hyp.rom.moderator.adj.sets%>%
  emmprep(at=list(avg.rir=0:23,
                  adj.sets.week=c(8,16,24)))%>%
  emmeans(~avg.rir+adj.sets.week,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.hyp.rom.moderator.adj.sets)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

primary.hyp.rom.moderator.adj.sets.contrast=emmprep(primary.hyp.rom.moderator.adj.sets,
                                                    at=list(avg.rir=primary.hyp.mean.rir,
                                                            adj.sets.week=primary.hyp.mean.adj.sets.week))%>%
  emmeans(~avg.rir+adj.sets.week,
          weights = "prop")%>%
  contrast("consec",interaction=TRUE)%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

primary.hyp.rom.moderator.adj.sets.contrast.label=paste0('Interaction Contrast (Adjusted Direct Sets) = ',
                                                         bquote(.(round(primary.hyp.rom.moderator.adj.sets.contrast[[3]],5))),
                                                         " [95% CI: ",bquote(.(round(primary.hyp.rom.moderator.adj.sets.contrast[[6]],5))),
                                                         ", ",bquote(.(round(primary.hyp.rom.moderator.adj.sets.contrast[[7]],5))),"]")

primary.hyp.rom.moderator.adj.sets.plot=ggplot(data = primary.data.hyp%>%
                                                 mutate(adj.sets.week=ifelse(adj.sets.week<12,"8 Sets",
                                                                             ifelse(adj.sets.week>12&adj.sets.week<20,"16 Sets","24 Sets"))),
                                               aes(x=avg.rir,
                                                   y=yi.rom.exp))+THEME+theme(plot.caption = element_text(face = "italic"))+
  facet_wrap(~fct_relevel(adj.sets.week,c("8 Sets","16 Sets","24 Sets")))+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Adjusted Direct Sets Per Week Per Muscle Group x Estimated RIR",
       caption = primary.hyp.rom.moderator.adj.sets.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = primary.hyp.rom.moderator.adj.sets.emms%>%
                mutate(adj.sets.week=ifelse(adj.sets.week==8,"8 Sets",
                                            ifelse(adj.sets.week==16,"16 Sets","24 Sets"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = primary.hyp.rom.moderator.adj.sets.emms%>%
                mutate(adj.sets.week=ifelse(adj.sets.week==8,"8 Sets",
                                            ifelse(adj.sets.week==16,"16 Sets","24 Sets"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = primary.hyp.rom.moderator.adj.sets.emms%>%
              mutate(adj.sets.week=ifelse(adj.sets.week==8,"8 Sets",
                                          ifelse(adj.sets.week==16,"16 Sets","24 Sets"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### rom frequency
primary.hyp.rom.moderator.frequency.per.muscle=primary.hyp.rom.best.model.prep%>%
  update(~. + avg.rir*frequency.per.muscle)%>%
  robust.rma.mv(cluster = primary.data.hyp$study)

primary.hyp.rom.moderator.frequency.per.muscle.emms=primary.hyp.rom.moderator.frequency.per.muscle%>%
  emmprep(at=list(avg.rir=0:23,
                  frequency.per.muscle=c(1,2,3)))%>%
  emmeans(~avg.rir+frequency.per.muscle,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.hyp.rom.moderator.frequency.per.muscle)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

primary.hyp.rom.moderator.frequency.per.muscle.contrast=emmprep(primary.hyp.rom.moderator.frequency.per.muscle,
                                                                at=list(avg.rir=primary.hyp.mean.rir,
                                                                        frequency.per.muscle=primary.hyp.mean.frequency.per.muscle))%>%
  emmeans(~avg.rir+frequency.per.muscle,
          weights = "prop")%>%
  contrast("consec",interaction=TRUE)%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

primary.hyp.rom.moderator.frequency.per.muscle.contrast.label=paste0('Interaction Contrast (Direct Frequency) = ',
                                                                     bquote(.(round(primary.hyp.rom.moderator.frequency.per.muscle.contrast[[3]],5))),
                                                                     " [95% CI: ",bquote(.(round(primary.hyp.rom.moderator.frequency.per.muscle.contrast[[6]],5))),
                                                                     ", ",bquote(.(round(primary.hyp.rom.moderator.frequency.per.muscle.contrast[[7]],5))),"]")

primary.hyp.rom.moderator.frequency.per.muscle.plot=ggplot(data = primary.data.hyp%>%
                                                             mutate(frequency.per.muscle=ifelse(frequency.per.muscle<1.5,"1x Per Week",
                                                                                                ifelse(frequency.per.muscle>1.5&frequency.per.muscle<2.5,"2x Per Week","3x Per Week"))),
                                                           aes(x=avg.rir,
                                                               y=yi.rom.exp))+THEME+theme(plot.caption = element_text(face = "italic"))+
  facet_wrap(~frequency.per.muscle)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Direct Frequency Per Week Per Muscle Group x Estimated RIR",
       caption = primary.hyp.rom.moderator.frequency.per.muscle.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = primary.hyp.rom.moderator.frequency.per.muscle.emms%>%
                mutate(frequency.per.muscle=ifelse(frequency.per.muscle==1,"1x Per Week",
                                                   ifelse(frequency.per.muscle==2,"2x Per Week","3x Per Week"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = primary.hyp.rom.moderator.frequency.per.muscle.emms%>%
                mutate(frequency.per.muscle=ifelse(frequency.per.muscle==1,"1x Per Week",
                                                   ifelse(frequency.per.muscle==2,"2x Per Week","3x Per Week"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = primary.hyp.rom.moderator.frequency.per.muscle.emms%>%
              mutate(frequency.per.muscle=ifelse(frequency.per.muscle==1,"1x Per Week",
                                                 ifelse(frequency.per.muscle==2,"2x Per Week","3x Per Week"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### rom age
primary.hyp.rom.moderator.age=primary.hyp.rom.best.model.prep%>%
  update(~. + avg.rir*age)%>%
  robust.rma.mv(cluster = primary.data.hyp$study)

primary.hyp.rom.moderator.age.emms=primary.hyp.rom.moderator.age%>%
  emmprep(at=list(avg.rir=0:23,
                  age=c(20,40,60)))%>%
  emmeans(~avg.rir+age,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.hyp.rom.moderator.age)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

primary.hyp.rom.moderator.age.contrast=emmprep(primary.hyp.rom.moderator.age,
                                               at=list(avg.rir=primary.hyp.mean.rir,
                                                       age=primary.hyp.mean.age))%>%
  emmeans(~avg.rir+age,
          weights = "prop")%>%
  contrast("consec",interaction=TRUE)%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

primary.hyp.rom.moderator.age.contrast.label=paste0('Interaction Contrast (Age) = ',
                                                    bquote(.(round(primary.hyp.rom.moderator.age.contrast[[3]],5))),
                                                    " [95% CI: ",bquote(.(round(primary.hyp.rom.moderator.age.contrast[[6]],5))),
                                                    ", ",bquote(.(round(primary.hyp.rom.moderator.age.contrast[[7]],5))),"]")

primary.hyp.rom.moderator.age.plot=ggplot(data = primary.data.hyp%>%
                                            mutate(age=ifelse(age<30,"20 Years Old",
                                                              ifelse(age>30&age<50,"40 Years Old","60 Years Old"))),
                                          aes(x=avg.rir,
                                              y=yi.rom.exp))+THEME+theme(plot.caption = element_text(face = "italic"))+
  facet_wrap(~age)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Age x Estimated RIR",
       caption = primary.hyp.rom.moderator.age.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = primary.hyp.rom.moderator.age.emms%>%
                mutate(age=ifelse(age==20,"20 Years Old",
                                  ifelse(age==40,"40 Years Old","60 Years Old"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = primary.hyp.rom.moderator.age.emms%>%
                mutate(age=ifelse(age==20,"20 Years Old",
                                  ifelse(age==40,"40 Years Old","60 Years Old"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = primary.hyp.rom.moderator.age.emms%>%
              mutate(age=ifelse(age==20,"20 Years Old",
                                ifelse(age==40,"40 Years Old","60 Years Old"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### rom sex
primary.hyp.rom.moderator.sex.percent.male=primary.hyp.rom.best.model.prep%>%
  update(~. + avg.rir*sex.percent.male)%>%
  robust.rma.mv(cluster = primary.data.hyp$study)

primary.hyp.rom.moderator.sex.percent.male.emms=primary.hyp.rom.moderator.sex.percent.male%>%
  emmprep(at=list(avg.rir=0:23,
                  sex.percent.male=c(0,100)))%>%
  emmeans(~avg.rir+sex.percent.male,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.hyp.rom.moderator.sex.percent.male)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

primary.hyp.rom.moderator.sex.percent.male.contrast=emmprep(primary.hyp.rom.moderator.sex.percent.male,
                                                            at=list(avg.rir=primary.hyp.mean.rir,
                                                                    sex.percent.male=primary.hyp.mean.sex))%>%
  emmeans(~avg.rir+sex.percent.male,
          weights = "prop")%>%
  contrast("consec",interaction=TRUE)%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

primary.hyp.rom.moderator.sex.percent.male.contrast.label=paste0('Interaction Contrast (Sex) = ',
                                                                 bquote(.(round(primary.hyp.rom.moderator.sex.percent.male.contrast[[3]],5))),
                                                                 " [95% CI: ",bquote(.(round(primary.hyp.rom.moderator.sex.percent.male.contrast[[6]],5))),
                                                                 ", ",bquote(.(round(primary.hyp.rom.moderator.sex.percent.male.contrast[[7]],5))),"]")

primary.hyp.rom.moderator.sex.percent.male.plot=ggplot(data = primary.data.hyp%>%
                                                         mutate(sex.percent.male=ifelse(sex.percent.male<50,"Female","Male")),
                                                       aes(x=avg.rir,
                                                           y=yi.rom.exp))+THEME+theme(plot.caption = element_text(face = "italic"))+
  facet_wrap(~sex.percent.male)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Sex (% Male) x Estimated RIR",
       caption = primary.hyp.rom.moderator.sex.percent.male.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = primary.hyp.rom.moderator.sex.percent.male.emms%>%
                mutate(sex.percent.male=ifelse(sex.percent.male==0,"Female","Male")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = primary.hyp.rom.moderator.sex.percent.male.emms%>%
                mutate(sex.percent.male=ifelse(sex.percent.male==0,"Female","Male")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = primary.hyp.rom.moderator.sex.percent.male.emms%>%
              mutate(sex.percent.male=ifelse(sex.percent.male==0,"Female","Male")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### rom study design
primary.hyp.rom.moderator.within.between.design=primary.hyp.rom.best.model.prep%>%
  update(~. + avg.rir*within.between.design)%>%
  robust.rma.mv(cluster = primary.data.hyp$study)

primary.hyp.rom.moderator.within.between.design.emms=primary.hyp.rom.moderator.within.between.design%>%emmprep(at=list(avg.rir=0:23))%>%
  emmeans(~avg.rir+within.between.design,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.hyp.rom.moderator.within.between.design)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

primary.hyp.rom.moderator.within.between.design.contrast=emmprep(primary.hyp.rom.moderator.within.between.design,
                                                                 at=list(avg.rir=primary.hyp.mean.rir))%>%
  emmeans(~avg.rir+within.between.design,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

primary.hyp.rom.moderator.within.between.design.contrast.label=paste0('Interaction Contrast (Between - Within) = ',
                                                                      bquote(.(round(primary.hyp.rom.moderator.within.between.design.contrast[[1,3]],5))),
                                                                      " [95% CI: ",bquote(.(round(primary.hyp.rom.moderator.within.between.design.contrast[[1,6]],5))),
                                                                      ", ",bquote(.(round(primary.hyp.rom.moderator.within.between.design.contrast[[1,7]],5))),"]")


primary.hyp.rom.moderator.within.between.design.plot=ggplot(data = primary.data.hyp%>%
                                                              mutate(within.between.design=ifelse(within.between.design=="Between","Between Subjects","Within Subjects")),
                                                            aes(x=avg.rir,
                                                                y=yi.rom.exp))+THEME+theme(legend.position = "right",
                                                                                           plot.caption = element_text(face = "italic"))+
  facet_wrap(~within.between.design)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Study Design x Estimated RIR",
       caption = primary.hyp.rom.moderator.within.between.design.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = primary.hyp.rom.moderator.within.between.design.emms%>%
                mutate(within.between.design=ifelse(within.between.design=="Between","Between Subjects","Within Subjects")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = primary.hyp.rom.moderator.within.between.design.emms%>%
                mutate(within.between.design=ifelse(within.between.design=="Between","Between Subjects","Within Subjects")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = primary.hyp.rom.moderator.within.between.design.emms%>%
              mutate(within.between.design=ifelse(within.between.design=="Between","Between Subjects","Within Subjects")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### rom upper / lower outcome
primary.hyp.rom.moderator.upper.lower.other=primary.hyp.rom.best.model.prep%>%
  update(~. + avg.rir*upper.lower.other)%>%
  robust.rma.mv(cluster = primary.data.hyp$study)

primary.hyp.rom.moderator.upper.lower.other.emms=primary.hyp.rom.moderator.upper.lower.other%>%emmprep(at=list(avg.rir=0:23))%>%
  emmeans(~avg.rir+upper.lower.other,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.hyp.rom.moderator.upper.lower.other)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

primary.hyp.rom.moderator.upper.lower.other.contrast=emmprep(primary.hyp.rom.moderator.upper.lower.other,
                                                             at=list(avg.rir=primary.hyp.mean.rir))%>%
  emmeans(~avg.rir+upper.lower.other,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

primary.hyp.rom.moderator.upper.lower.other.contrast.label=paste0('Interaction Contrast (Lower - Upper) = ',
                                                                  bquote(.(round(primary.hyp.rom.moderator.upper.lower.other.contrast[[1,3]],5))),
                                                                  " [95% CI: ",bquote(.(round(primary.hyp.rom.moderator.upper.lower.other.contrast[[1,6]],5))),
                                                                  ", ",bquote(.(round(primary.hyp.rom.moderator.upper.lower.other.contrast[[1,7]],5))),"]")


primary.hyp.rom.moderator.upper.lower.other.plot=ggplot(data = primary.data.hyp%>%
                                                          mutate(upper.lower.other=ifelse(upper.lower.other=="Lower","Lower Body Outcomes","Upper Body Outcomes")),
                                                        aes(x=avg.rir,
                                                            y=yi.rom.exp))+THEME+theme(legend.position = "right",
                                                                                       plot.caption = element_text(face = "italic"))+
  facet_wrap(~upper.lower.other)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Upper/Lower Body Outcome x Estimated RIR",
       caption = primary.hyp.rom.moderator.upper.lower.other.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = primary.hyp.rom.moderator.upper.lower.other.emms%>%
                mutate(upper.lower.other=ifelse(upper.lower.other=="Lower","Lower Body Outcomes","Upper Body Outcomes")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = primary.hyp.rom.moderator.upper.lower.other.emms%>%
                mutate(upper.lower.other=ifelse(upper.lower.other=="Lower","Lower Body Outcomes","Upper Body Outcomes")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = primary.hyp.rom.moderator.upper.lower.other.emms%>%
              mutate(upper.lower.other=ifelse(upper.lower.other=="Lower","Lower Body Outcomes","Upper Body Outcomes")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### rom exercise selection
primary.hyp.rom.moderator.train.exercise=primary.hyp.rom.best.model.prep%>%
  update(~. + avg.rir*train.exercise)%>%
  robust.rma.mv(cluster = primary.data.hyp$study)

primary.hyp.rom.moderator.train.exercise.emms=primary.hyp.rom.moderator.train.exercise%>%emmprep(at=list(avg.rir=0:23))%>%
  emmeans(~avg.rir+train.exercise,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.hyp.rom.moderator.train.exercise)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

primary.hyp.rom.moderator.train.exercise.contrast=emmprep(primary.hyp.rom.moderator.train.exercise,
                                                          at=list(avg.rir=primary.hyp.mean.rir))%>%
  emmeans(~avg.rir+train.exercise,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

primary.hyp.rom.moderator.train.exercise.contrast.label=paste0('Interaction Contrast (Both - Multi-joint) = ',
                                                               bquote(.(round(primary.hyp.rom.moderator.train.exercise.contrast[[1,3]],5))),
                                                               " [95% CI: ",bquote(.(round(primary.hyp.rom.moderator.train.exercise.contrast[[1,6]],5))),
                                                               ", ",bquote(.(round(primary.hyp.rom.moderator.train.exercise.contrast[[1,7]],5))),"]\n",
                                                               'Interaction Contrast (Both - Single-joint) = ',
                                                               bquote(.(round(primary.hyp.rom.moderator.train.exercise.contrast[[2,3]],5))),
                                                               " [95% CI: ",bquote(.(round(primary.hyp.rom.moderator.train.exercise.contrast[[2,6]],5))),
                                                               ", ",bquote(.(round(primary.hyp.rom.moderator.train.exercise.contrast[[2,7]],5))),"]\n",
                                                               'Interaction Contrast (Multi-joint - Single-joint) = ',
                                                               bquote(.(round(primary.hyp.rom.moderator.train.exercise.contrast[[3,3]],5))),
                                                               " [95% CI: ",bquote(.(round(primary.hyp.rom.moderator.train.exercise.contrast[[3,6]],5))),
                                                               ", ",bquote(.(round(primary.hyp.rom.moderator.train.exercise.contrast[[3,7]],5))),"]")


primary.hyp.rom.moderator.train.exercise.plot=ggplot(data = primary.data.hyp%>%
                                                       mutate(train.exercise=ifelse(train.exercise=="Both","Both",
                                                                                    ifelse(train.exercise=="Multi-joint","Multi-joint Exercises","Single-joint Exercises"))),
                                                     aes(x=avg.rir,
                                                         y=yi.rom.exp))+THEME+theme(legend.position = "right",
                                                                                    plot.caption = element_text(face = "italic"))+
  facet_wrap(~train.exercise)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Exercise Selection x Estimated RIR",
       caption = primary.hyp.rom.moderator.train.exercise.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = primary.hyp.rom.moderator.train.exercise.emms%>%
                mutate(train.exercise=ifelse(train.exercise=="Both","Both",
                                             ifelse(train.exercise=="Multi-joint","Multi-joint Exercises","Single-joint Exercises"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = primary.hyp.rom.moderator.train.exercise.emms%>%
                mutate(train.exercise=ifelse(train.exercise=="Both","Both",
                                             ifelse(train.exercise=="Multi-joint","Multi-joint Exercises","Single-joint Exercises"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = primary.hyp.rom.moderator.train.exercise.emms%>%
              mutate(train.exercise=ifelse(train.exercise=="Both","Both",
                                           ifelse(train.exercise=="Multi-joint","Multi-joint Exercises","Single-joint Exercises"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### rom concurrent training
primary.hyp.rom.moderator.formal.cardio.intervention=primary.hyp.rom.best.model.prep%>%
  update(~. + avg.rir*formal.cardio.intervention)%>%
  robust.rma.mv(cluster = primary.data.hyp$study)

primary.hyp.rom.moderator.formal.cardio.intervention.emms=primary.hyp.rom.moderator.formal.cardio.intervention%>%emmprep(at=list(avg.rir=0:23))%>%
  emmeans(~avg.rir+formal.cardio.intervention,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.hyp.rom.moderator.formal.cardio.intervention)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

primary.hyp.rom.moderator.formal.cardio.intervention.contrast=emmprep(primary.hyp.rom.moderator.formal.cardio.intervention,
                                                                      at=list(avg.rir=primary.hyp.mean.rir))%>%
  emmeans(~avg.rir+formal.cardio.intervention,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

primary.hyp.rom.moderator.formal.cardio.intervention.contrast.label=paste0('Interaction Contrast (RT Only - Concurrent) = ',
                                                                           bquote(.(round(primary.hyp.rom.moderator.formal.cardio.intervention.contrast[[1,3]],5))),
                                                                           " [95% CI: ",bquote(.(round(primary.hyp.rom.moderator.formal.cardio.intervention.contrast[[1,6]],5))),
                                                                           ", ",bquote(.(round(primary.hyp.rom.moderator.formal.cardio.intervention.contrast[[1,7]],5))),"]")


primary.hyp.rom.moderator.formal.cardio.intervention.plot=ggplot(data = primary.data.hyp%>%
                                                                   mutate(formal.cardio.intervention=ifelse(formal.cardio.intervention=="No","Only Resistance Training","Concurrent Training")),
                                                                 aes(x=avg.rir,
                                                                     y=yi.rom.exp))+THEME+theme(legend.position = "right",
                                                                                                plot.caption = element_text(face = "italic"))+
  facet_wrap(~formal.cardio.intervention)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Concurrent Training x Estimated RIR",
       caption = primary.hyp.rom.moderator.formal.cardio.intervention.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = primary.hyp.rom.moderator.formal.cardio.intervention.emms%>%
                mutate(formal.cardio.intervention=ifelse(formal.cardio.intervention=="No","Only Resistance Training","Concurrent Training")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = primary.hyp.rom.moderator.formal.cardio.intervention.emms%>%
                mutate(formal.cardio.intervention=ifelse(formal.cardio.intervention=="No","Only Resistance Training","Concurrent Training")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = primary.hyp.rom.moderator.formal.cardio.intervention.emms%>%
              mutate(formal.cardio.intervention=ifelse(formal.cardio.intervention=="No","Only Resistance Training","Concurrent Training")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### rom progressive overload
primary.hyp.rom.moderator.progression=primary.hyp.rom.best.model.prep%>%
  update(~. + avg.rir*progression)%>%
  robust.rma.mv(cluster = primary.data.hyp$study)

primary.hyp.rom.moderator.progression.emms=primary.hyp.rom.moderator.progression%>%emmprep(at=list(avg.rir=0:23))%>%
  emmeans(~avg.rir+progression,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.hyp.rom.moderator.progression)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

primary.hyp.rom.moderator.progression.contrast=emmprep(primary.hyp.rom.moderator.progression,
                                                       at=list(avg.rir=primary.hyp.mean.rir))%>%
  emmeans(~avg.rir+progression,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

primary.hyp.rom.moderator.progression.contrast.label=paste0('Interaction Contrast (No Overload - Progressive Overload) = ',
                                                            bquote(.(round(primary.hyp.rom.moderator.progression.contrast[[1,3]],5))),
                                                            " [95% CI: ",bquote(.(round(primary.hyp.rom.moderator.progression.contrast[[1,6]],5))),
                                                            ", ",bquote(.(round(primary.hyp.rom.moderator.progression.contrast[[1,7]],5))),"]")


primary.hyp.rom.moderator.progression.plot=ggplot(data = primary.data.hyp%>%
                                                    mutate(progression=ifelse(progression=="No","No Progressive Overload","Progressive Overload")),
                                                  aes(x=avg.rir,
                                                      y=yi.rom.exp))+THEME+theme(legend.position = "right",
                                                                                 plot.caption = element_text(face = "italic"))+
  facet_wrap(~progression)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Progressive Overload x Estimated RIR",
       caption = primary.hyp.rom.moderator.progression.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = primary.hyp.rom.moderator.progression.emms%>%
                mutate(progression=ifelse(progression=="No","No Progressive Overload","Progressive Overload")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = primary.hyp.rom.moderator.progression.emms%>%
                mutate(progression=ifelse(progression=="No","No Progressive Overload","Progressive Overload")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = primary.hyp.rom.moderator.progression.emms%>%
              mutate(progression=ifelse(progression=="No","No Progressive Overload","Progressive Overload")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### rom failure definition
primary.hyp.rom.moderator.failure.definition=primary.hyp.rom.best.model.prep%>%
  update(~. + avg.rir*failure.definition)%>%
  robust.rma.mv(cluster = primary.data.hyp$study)

primary.hyp.rom.moderator.failure.definition.emms=primary.hyp.rom.moderator.failure.definition%>%emmprep(at=list(avg.rir=0:23))%>%
  emmeans(~avg.rir+failure.definition,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.hyp.rom.moderator.failure.definition)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

primary.hyp.rom.moderator.failure.definition.contrast=emmprep(primary.hyp.rom.moderator.failure.definition,
                                                              at=list(avg.rir=primary.hyp.mean.rir))%>%
  emmeans(~avg.rir+failure.definition,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

primary.hyp.rom.moderator.failure.definition.contrast.label=paste0('Interaction Contrast (No Definition - Definition Provided) = ',
                                                                   bquote(.(round(primary.hyp.rom.moderator.failure.definition.contrast[[1,3]],5))),
                                                                   " [95% CI: ",bquote(.(round(primary.hyp.rom.moderator.failure.definition.contrast[[1,6]],5))),
                                                                   ", ",bquote(.(round(primary.hyp.rom.moderator.failure.definition.contrast[[1,7]],5))),"]")


primary.hyp.rom.moderator.failure.definition.plot=ggplot(data = primary.data.hyp%>%
                                                           mutate(failure.definition=ifelse(failure.definition=="No","No Failure Definition","Failure Definition Provided")),
                                                         aes(x=avg.rir,
                                                             y=yi.rom.exp))+THEME+theme(legend.position = "right",
                                                                                        plot.caption = element_text(face = "italic"))+
  facet_wrap(~failure.definition)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Failure Definition x Estimated RIR",
       caption = primary.hyp.rom.moderator.failure.definition.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = primary.hyp.rom.moderator.failure.definition.emms%>%
                mutate(failure.definition=ifelse(failure.definition=="No","No Failure Definition","Failure Definition Provided")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = primary.hyp.rom.moderator.failure.definition.emms%>%
                mutate(failure.definition=ifelse(failure.definition=="No","No Failure Definition","Failure Definition Provided")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = primary.hyp.rom.moderator.failure.definition.emms%>%
              mutate(failure.definition=ifelse(failure.definition=="No","No Failure Definition","Failure Definition Provided")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### rom set structure
primary.hyp.rom.moderator.alternative.set.structure=primary.hyp.rom.best.model.prep%>%
  update(~. + avg.rir*alternative.set.structure)%>%
  robust.rma.mv(cluster = primary.data.hyp$study)

primary.hyp.rom.moderator.alternative.set.structure.emms=primary.hyp.rom.moderator.alternative.set.structure%>%emmprep(at=list(avg.rir=0:23))%>%
  emmeans(~avg.rir+alternative.set.structure,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.hyp.rom.moderator.alternative.set.structure)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

primary.hyp.rom.moderator.alternative.set.structure.contrast=emmprep(primary.hyp.rom.moderator.alternative.set.structure,
                                                                     at=list(avg.rir=primary.hyp.mean.rir))%>%
  emmeans(~avg.rir+alternative.set.structure,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

primary.hyp.rom.moderator.alternative.set.structure.contrast.label=paste0('Interaction Contrast (Traditional - Alternative) = ',
                                                                          bquote(.(round(primary.hyp.rom.moderator.alternative.set.structure.contrast[[1,3]],5))),
                                                                          " [95% CI: ",bquote(.(round(primary.hyp.rom.moderator.alternative.set.structure.contrast[[1,6]],5))),
                                                                          ", ",bquote(.(round(primary.hyp.rom.moderator.alternative.set.structure.contrast[[1,7]],5))),"]")


primary.hyp.rom.moderator.alternative.set.structure.plot=ggplot(data = primary.data.hyp%>%
                                                                  mutate(alternative.set.structure=ifelse(alternative.set.structure=="No","Traditional Set Structure","Alternative Set Structure")),
                                                                aes(x=avg.rir,
                                                                    y=yi.rom.exp))+THEME+theme(legend.position = "right",
                                                                                               plot.caption = element_text(face = "italic"))+
  facet_wrap(~alternative.set.structure)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Set Structure x Estimated RIR",
       caption = primary.hyp.rom.moderator.alternative.set.structure.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = primary.hyp.rom.moderator.alternative.set.structure.emms%>%
                mutate(alternative.set.structure=ifelse(alternative.set.structure=="No","Traditional Set Structure","Alternative Set Structure")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = primary.hyp.rom.moderator.alternative.set.structure.emms%>%
                mutate(alternative.set.structure=ifelse(alternative.set.structure=="No","Traditional Set Structure","Alternative Set Structure")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = primary.hyp.rom.moderator.alternative.set.structure.emms%>%
              mutate(alternative.set.structure=ifelse(alternative.set.structure=="No","Traditional Set Structure","Alternative Set Structure")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### rom intended velocity
primary.hyp.rom.moderator.train.intent=primary.hyp.rom.best.model.prep%>%
  update(~. + avg.rir*train.intent)%>%
  robust.rma.mv(cluster = primary.data.hyp$study)

primary.hyp.rom.moderator.train.intent.emms=primary.hyp.rom.moderator.train.intent%>%emmprep(at=list(avg.rir=0:23))%>%
  emmeans(~avg.rir+train.intent,
          weights = "prop")%>%
  pred_interval_esmeans(model=primary.hyp.rom.moderator.train.intent)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

primary.hyp.rom.moderator.train.intent.contrast=emmprep(primary.hyp.rom.moderator.train.intent,
                                                        at=list(avg.rir=primary.hyp.mean.rir))%>%
  emmeans(~avg.rir+train.intent,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

primary.hyp.rom.moderator.train.intent.contrast.label=paste0('Interaction Contrast (Submaximal - Maximal) = ',
                                                             bquote(.(round(primary.hyp.rom.moderator.train.intent.contrast[[1,3]],5))),
                                                             " [95% CI: ",bquote(.(round(primary.hyp.rom.moderator.train.intent.contrast[[1,6]],5))),
                                                             ", ",bquote(.(round(primary.hyp.rom.moderator.train.intent.contrast[[1,7]],5))),"]")


primary.hyp.rom.moderator.train.intent.plot=ggplot(data = primary.data.hyp%>%
                                                     mutate(train.intent=ifelse(train.intent=="No","Submaximal Intended Concentric Velocity","Maximal Intended Concentric Velocity")),
                                                   aes(x=avg.rir,
                                                       y=yi.rom.exp))+THEME+theme(legend.position = "right",
                                                                                  plot.caption = element_text(face = "italic"))+
  facet_wrap(~train.intent)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Intended Concentric Velocity x Estimated RIR",
       caption = primary.hyp.rom.moderator.train.intent.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = primary.hyp.rom.moderator.train.intent.emms%>%
                mutate(train.intent=ifelse(train.intent=="No","Submaximal Intended Concentric Velocity","Maximal Intended Concentric Velocity")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = primary.hyp.rom.moderator.train.intent.emms%>%
                mutate(train.intent=ifelse(train.intent=="No","Submaximal Intended Concentric Velocity","Maximal Intended Concentric Velocity")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = primary.hyp.rom.moderator.train.intent.emms%>%
              mutate(train.intent=ifelse(train.intent=="No","Submaximal Intended Concentric Velocity","Maximal Intended Concentric Velocity")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### combine rom plots
primary.hyp.rom.moderator.doc=primary.hyp.rom.moderator.load.plot/plot_spacer()/primary.hyp.rom.moderator.volume.method.plot/plot_spacer()/
  primary.hyp.rom.moderator.train.status.plot/plot_spacer()/primary.hyp.rom.moderator.weeks.plot/plot_spacer()/
  primary.hyp.rom.moderator.adj.sets.plot/plot_spacer()/primary.hyp.rom.moderator.frequency.per.muscle.plot/plot_spacer()/
  primary.hyp.rom.moderator.age.plot/plot_spacer()/primary.hyp.rom.moderator.sex.percent.male.plot/plot_spacer()/
  primary.hyp.rom.moderator.within.between.design.plot/plot_spacer()/primary.hyp.rom.moderator.upper.lower.other.plot/plot_spacer()/
  primary.hyp.rom.moderator.train.exercise.plot/plot_spacer()/primary.hyp.rom.moderator.formal.cardio.intervention.plot/plot_spacer()/
  primary.hyp.rom.moderator.progression.plot/plot_spacer()/primary.hyp.rom.moderator.failure.definition.plot/plot_spacer()/
  primary.hyp.rom.moderator.alternative.set.structure.plot/plot_spacer()/primary.hyp.rom.moderator.train.intent.plot&
  plot_annotation(title = "Moderator Analyses from Linear Multilevel Model for Muscle Hypertrophy (Response Ratio)")&
  THEME

ggsave(primary.hyp.rom.moderator.doc,device = "pdf",filename="primary.hyp.rom.moderator.doc.pdf",width = 12,height = 136.5,limitsize = FALSE,
       path = "/Users/zacrobinson/Documents/For Publication PTF Meta/Plots")

### SECONDARY ANALYSIS (<10 RIR) --------------------------------------------
# CHARACTERISTICS ---------------------------------------------------------
#Visual Descriptions of Training Interventions
secondary.densitites.volume<-ggplot(data = secondary.data,aes(x=adj.sets.week,color=outcome,fill=outcome))+THEME+
  theme(axis.text.y = element_text(face = "bold",angle = 90,size = 8,hjust = 0.5))+
  scale_x_continuous(breaks = c(0,2,4,6,8,10))+
  scale_color_manual(values = c("#0000a7","#c1272d"))+
  scale_fill_manual(values = c("#0000a7","#c1272d"))+
  geom_density_ridges(aes(y=outcome),
                      jittered_points = TRUE,
                      position = "raincloud",
                      alpha = 0.7, scale = 0.5)+
  labs(x="Direct Adjusted Sets Per Week Per Muscle Group",y="Density")

secondary.densitites.load<-ggplot(data = secondary.data,aes(x=load.set,color=outcome,fill=outcome))+THEME+
  theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank())+
  scale_x_continuous(breaks = c(30,40,50,60,70,80,90))+
  scale_color_manual(values = c("#0000a7","#c1272d"))+
  scale_fill_manual(values = c("#0000a7","#c1272d"))+
  geom_density_ridges(aes(y=outcome),
                      jittered_points = TRUE,
                      position = "raincloud",
                      alpha = 0.7, scale = 0.5)+
  labs(x="Load (% of 1RM) Per Set",y="Density")

secondary.densitites.frequency<-ggplot(data = secondary.data,aes(x=frequency.per.muscle,color=outcome,fill=outcome))+THEME+
  theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank())+
  scale_x_continuous(breaks = c(1,2,3,4,5,6))+
  scale_color_manual(values = c("#0000a7","#c1272d"))+
  scale_fill_manual(values = c("#0000a7","#c1272d"))+
  geom_density_ridges(aes(y=outcome),
                      jittered_points = TRUE,
                      position = "raincloud",
                      alpha = 0.7, scale = 0.5)+
  labs(x="Direct Frequency (Sessions Per Week Per Muscle Group)",y="Density")

secondary.density.plot.final<-secondary.densitites.volume+secondary.densitites.load+secondary.densitites.frequency&
  plot_annotation(title = "Visual Summaries of Training Interventions Included in Meta-regression Models <10 RIR")&
  theme(plot.title = element_text(size = 15,hjust = 0.5))&THEME

ggsave(secondary.density.plot.final,device = "pdf",filename="secondary.density.plot.final.pdf",width = 12,height = 6,
       path = "/Users/zacrobinson/Documents/For Publication PTF Meta/Plots")

#Study Characteristics
secondary.chars<-c('train.status','nutrition.controlled','formal.cardio.intervention','progression','failure.definition',
                 'alternative.set.structure','set.rep.equated','upper.lower.other','train.exercise','within.between.design','train.intent')

for (x in secondary.chars) {
  list<-secondary.data%>%
    group_by(outcome,.data[[x]])%>%
    count(.data[[x]])%>%
    spread(outcome,n)%>%
    mutate(Variable=x)%>%
    rename(Level=.data[[x]])%>%
    relocate(Variable)
  
  assign((paste0("secondary.char.table.",x)),list)
}

secondary.char.table.combined<-rbind(secondary.char.table.train.status,
                                   secondary.char.table.nutrition.controlled,
                                   secondary.char.table.formal.cardio.intervention,
                                   secondary.char.table.progression,
                                   secondary.char.table.failure.definition,
                                   secondary.char.table.alternative.set.structure,
                                   secondary.char.table.set.rep.equated,
                                   secondary.char.table.upper.lower.other,
                                   secondary.char.table.train.exercise,
                                   secondary.char.table.within.between.design,
                                   secondary.char.table.train.intent)

secondary.char.table.combined$Variable=c(rep("Training Status",2),
                                       rep("Nutrition Control",2),
                                       rep("Concurrent Training Intervention",2),
                                       rep("Employed Progressive Overload",2),
                                       rep("Failure Definition Provided",2),
                                       rep("Alternative Set Structures",2),
                                       rep("Method of Equating Volume",3),
                                       rep("Upper or Lower Body Outcome",2),
                                       rep("Training Exercise Selection",3),
                                       rep("Within or Between Participant Design",2),
                                       rep("Utilized Maximal Intended Concentric Velocity",2))

secondary.char.table.combined$Level=c("Trained","Untrained",
                                    "No","Yes",
                                    "No","Yes",
                                    "No","Yes",
                                    "No","Yes",
                                    "No","Yes",
                                    "Set","Repetition","Both",
                                    "Lower","Upper","Both",
                                    "Multi-joint","Single-joint",
                                    "Between","Within",
                                    "No","Yes")

secondary.char.table.combined.final<-secondary.char.table.combined%>%
  gt(groupname_col =  "Variable")%>%
  tab_header("Characteristics of Effects Included in Meta-regression Models with <10 RIR")%>%
  opt_stylize(style = 1,color = "gray")%>%
  tab_options(row_group.as_column = TRUE)%>%
  tab_style(style = cell_text(style = "italic"),locations = cells_row_groups())%>%
  tab_style(style = cell_text(align = "center"),locations = cells_column_labels())%>%
  tab_style(style = cell_text(weight  = "bold"),locations = cells_title())%>%
  tab_style(style = cell_text(align = "center"),locations = cells_body())

gtsave(secondary.char.table.combined.final,filename = "secondary.char.table.combined.final.png",
       path = "/Users/zacrobinson/Documents/For Publication PTF Meta/Tables")

# FIT AND COMPARE MODELS --------------------------------------------------------------

### strength
#### linear
secondary.str.g.linear.model.ri=rma.mv(yi,vi,data = secondary.data.str,
                                     mods = ~ avg.rir+load.set+set.rep.equated+weeks+train.status,
                                     random = list(~1|study,~1|group,~1|obs),
                                     method = "REML", test = "t",dfs = "contain")

secondary.str.g.linear.model.rs=rma.mv(yi,vi,data = secondary.data.str,
                                     mods = ~ avg.rir+load.set+set.rep.equated+weeks+train.status,
                                     random = list(~avg.rir|study,~1|group,~1|obs),
                                     method = "REML", test = "t",dfs = "contain",struct = "GEN")

secondary.str.g.linear.model.rs2=rma.mv(yi,vi,data = secondary.data.str,
                                      mods = ~ avg.rir+load.set+set.rep.equated+weeks+train.status,
                                      random = list(~avg.rir|study,~avg.rir|group,~1|obs),
                                      method = "REML", test = "t",dfs = "contain",struct = c("GEN","GEN"))

secondary.str.rom.linear.model.ri=rma.mv(yi.rom,vi.rom,data = secondary.data.str,
                                       mods = ~ avg.rir+load.set+set.rep.equated+weeks+train.status,
                                       random = list(~1|study,~1|group,~1|obs),
                                       method = "REML", test = "t",dfs = "contain")

secondary.str.rom.linear.model.rs=rma.mv(yi.rom,vi.rom,data = secondary.data.str,
                                       mods = ~ avg.rir+load.set+set.rep.equated+weeks+train.status,
                                       random = list(~avg.rir|study,~1|group,~1|obs),
                                       method = "REML", test = "t",dfs = "contain",struct = "GEN")

secondary.str.rom.linear.model.rs2=rma.mv(yi.rom,vi.rom,data = secondary.data.str,
                                        mods = ~ avg.rir+load.set+set.rep.equated+weeks+train.status,
                                        random = list(~avg.rir|study,~avg.rir|group,~1|obs),
                                        method = "REML", test = "t",dfs = "contain",struct = c("GEN","GEN"))

#### linear-log
secondary.str.g.log.model.ri=rma.mv(yi,vi,data = secondary.data.str,
                                  mods = ~ log1p(avg.rir)+load.set+set.rep.equated+weeks+train.status,
                                  random = list(~1|study,~1|group,~1|obs),
                                  method = "REML", test = "t",dfs = "contain")

secondary.str.g.log.model.rs=rma.mv(yi,vi,data = secondary.data.str,
                                  mods = ~ log1p(avg.rir)+load.set+set.rep.equated+weeks+train.status,
                                  random = list(~avg.rir|study,~1|group,~1|obs),
                                  method = "REML", test = "t",dfs = "contain",struct = "GEN")

secondary.str.g.log.model.rs2=rma.mv(yi,vi,data = secondary.data.str,
                                   mods = ~ log1p(avg.rir)+load.set+set.rep.equated+weeks+train.status,
                                   random = list(~avg.rir|study,~avg.rir|group,~1|obs),
                                   method = "REML", test = "t",dfs = "contain",struct = c("GEN","GEN"))

secondary.str.rom.log.model.ri=rma.mv(yi.rom,vi.rom,data = secondary.data.str,
                                    mods = ~ log1p(avg.rir)+load.set+set.rep.equated+weeks+train.status,
                                    random = list(~1|study,~1|group,~1|obs),
                                    method = "REML", test = "t",dfs = "contain")

secondary.str.rom.log.model.rs=rma.mv(yi.rom,vi.rom,data = secondary.data.str,
                                    mods = ~ log1p(avg.rir)+load.set+set.rep.equated+weeks+train.status,
                                    random = list(~avg.rir|study,~1|group,~1|obs),
                                    method = "REML", test = "t",dfs = "contain",struct = "GEN")

secondary.str.rom.log.model.rs2=rma.mv(yi.rom,vi.rom,data = secondary.data.str,
                                     mods = ~ log1p(avg.rir)+load.set+set.rep.equated+weeks+train.status,
                                     random = list(~avg.rir|study,~avg.rir|group,~1|obs),
                                     method = "REML", test = "t",dfs = "contain",struct = c("GEN","GEN"))

#### quadratic
secondary.str.g.quadratic.model.ri=rma.mv(yi,vi,data = secondary.data.str,
                                        mods = ~ poly(avg.rir,2,raw = TRUE)+load.set+set.rep.equated+weeks+train.status,
                                        random = list(~1|study,~1|group,~1|obs),
                                        method = "REML", test = "t",dfs = "contain")

secondary.str.g.quadratic.model.rs=rma.mv(yi,vi,data = secondary.data.str,
                                        mods = ~ poly(avg.rir,2,raw = TRUE)+load.set+set.rep.equated+weeks+train.status,
                                        random = list(~avg.rir|study,~1|group,~1|obs),
                                        method = "REML", test = "t",dfs = "contain",struct = "GEN")

secondary.str.g.quadratic.model.rs2=rma.mv(yi,vi,data = secondary.data.str,
                                         mods = ~ poly(avg.rir,2,raw = TRUE)+load.set+set.rep.equated+weeks+train.status,
                                         random = list(~avg.rir|study,~avg.rir|group,~1|obs),
                                         method = "REML", test = "t",dfs = "contain",struct = c("GEN","GEN"))

secondary.str.rom.quadratic.model.ri=rma.mv(yi.rom,vi.rom,data = secondary.data.str,
                                          mods = ~ poly(avg.rir,2,raw = TRUE)+load.set+set.rep.equated+weeks+train.status,
                                          random = list(~1|study,~1|group,~1|obs),
                                          method = "REML", test = "t",dfs = "contain")

secondary.str.rom.quadratic.model.rs=rma.mv(yi.rom,vi.rom,data = secondary.data.str,
                                          mods = ~ poly(avg.rir,2,raw = TRUE)+load.set+set.rep.equated+weeks+train.status,
                                          random = list(~avg.rir|study,~1|group,~1|obs),
                                          method = "REML", test = "t",dfs = "contain",struct = "GEN")

secondary.str.rom.quadratic.model.rs2=rma.mv(yi.rom,vi.rom,data = secondary.data.str,
                                           mods = ~ poly(avg.rir,2,raw = TRUE)+load.set+set.rep.equated+weeks+train.status,
                                           random = list(~avg.rir|study,~avg.rir|group,~1|obs),
                                           method = "REML", test = "t",dfs = "contain",struct = c("GEN","GEN"))

#### cubic
secondary.str.g.cubic.model.ri=rma.mv(yi,vi,data = secondary.data.str,
                                    mods = ~ poly(avg.rir,3,raw = TRUE)+load.set+set.rep.equated+weeks+train.status,
                                    random = list(~1|study,~1|group,~1|obs),
                                    method = "REML", test = "t",dfs = "contain")

secondary.str.g.cubic.model.rs=rma.mv(yi,vi,data = secondary.data.str,
                                    mods = ~ poly(avg.rir,3,raw = TRUE)+load.set+set.rep.equated+weeks+train.status,
                                    random = list(~avg.rir|study,~1|group,~1|obs),
                                    method = "REML", test = "t",dfs = "contain",struct = "GEN")

secondary.str.g.cubic.model.rs2=rma.mv(yi,vi,data = secondary.data.str,
                                     mods = ~ poly(avg.rir,3,raw = TRUE)+load.set+set.rep.equated+weeks+train.status,
                                     random = list(~avg.rir|study,~avg.rir|group,~1|obs),
                                     method = "REML", test = "t",dfs = "contain",struct = c("GEN","GEN"))

secondary.str.rom.cubic.model.ri=rma.mv(yi.rom,vi.rom,data = secondary.data.str,
                                      mods = ~ poly(avg.rir,3,raw = TRUE)+load.set+set.rep.equated+weeks+train.status,
                                      random = list(~1|study,~1|group,~1|obs),
                                      method = "REML", test = "t",dfs = "contain")

secondary.str.rom.cubic.model.rs=rma.mv(yi.rom,vi.rom,data = secondary.data.str,
                                      mods = ~ poly(avg.rir,3,raw = TRUE)+load.set+set.rep.equated+weeks+train.status,
                                      random = list(~avg.rir|study,~1|group,~1|obs),
                                      method = "REML", test = "t",dfs = "contain",struct = "GEN")

secondary.str.rom.cubic.model.rs2=rma.mv(yi.rom,vi.rom,data = secondary.data.str,
                                       mods = ~ poly(avg.rir,3,raw = TRUE)+load.set+set.rep.equated+weeks+train.status,
                                       random = list(~avg.rir|study,~avg.rir|group,~1|obs),
                                       method = "REML", test = "t",dfs = "contain",struct = c("GEN","GEN"))

#### linear-spline
secondary.str.g.spline.model.ri=rma.mv(yi,vi,data = secondary.data.str,
                                     mods = ~ bs(spline.rir,degree = 1,knots = 0)+load.set+set.rep.equated+weeks+train.status,
                                     random = list(~1|study,~1|group,~1|obs),
                                     method = "REML", test = "t",dfs = "contain")

secondary.str.g.spline.model.rs=rma.mv(yi,vi,data = secondary.data.str,
                                     mods = ~ bs(spline.rir,degree = 1,knots = 0)+load.set+set.rep.equated+weeks+train.status,
                                     random = list(~spline.rir|study,~1|group,~1|obs),
                                     method = "REML", test = "t",dfs = "contain",struct = "GEN")

secondary.str.g.spline.model.rs2=rma.mv(yi,vi,data = secondary.data.str,
                                      mods = ~ bs(spline.rir,degree = 1,knots = 0)+load.set+set.rep.equated+weeks+train.status,
                                      random = list(~spline.rir|study,~spline.rir|group,~1|obs),
                                      method = "REML", test = "t",dfs = "contain",struct = c("GEN","GEN"))

secondary.str.rom.spline.model.ri=rma.mv(yi.rom,vi.rom,data = secondary.data.str,
                                       mods = ~ bs(spline.rir,degree = 1,knots = 0)+load.set+set.rep.equated+weeks+train.status,
                                       random = list(~1|study,~1|group,~1|obs),
                                       method = "REML", test = "t",dfs = "contain")

secondary.str.rom.spline.model.rs=rma.mv(yi.rom,vi.rom,data = secondary.data.str,
                                       mods = ~ bs(spline.rir,degree = 1,knots = 0)+load.set+set.rep.equated+weeks+train.status,
                                       random = list(~spline.rir|study,~1|group,~1|obs),
                                       method = "REML", test = "t",dfs = "contain",struct = "GEN")

secondary.str.rom.spline.model.rs2=rma.mv(yi.rom,vi.rom,data = secondary.data.str,
                                        mods = ~ bs(spline.rir,degree = 1,knots = 0)+load.set+set.rep.equated+weeks+train.status,
                                        random = list(~spline.rir|study,~spline.rir|group,~1|obs),
                                        method = "REML", test = "t",dfs = "contain",struct = c("GEN","GEN"))

#### cubic-spline
secondary.str.g.rcs.model.ri=rma.mv(yi,vi,data = secondary.data.str,
                                  mods = ~ ns(avg.rir,knots = secondary.str.k,Boundary.knots=secondary.str.b)+load.set+set.rep.equated+weeks+train.status,
                                  random = list(~1|study,~1|group,~1|obs),
                                  method = "REML", test = "t",dfs = "contain")

secondary.str.g.rcs.model.rs=rma.mv(yi,vi,data = secondary.data.str,
                                  mods = ~ ns(avg.rir,knots = secondary.str.k,Boundary.knots = secondary.str.b)+load.set+set.rep.equated+weeks+train.status,
                                  random = list(~avg.rir|study,~1|group,~1|obs),
                                  method = "REML", test = "t",dfs = "contain",struct = "GEN")

secondary.str.g.rcs.model.rs2=rma.mv(yi,vi,data = secondary.data.str,
                                   mods = ~ ns(avg.rir,knots = secondary.str.k,Boundary.knots = secondary.str.b)+load.set+set.rep.equated+weeks+train.status,
                                   random = list(~avg.rir|study,~avg.rir|group,~1|obs),
                                   method = "REML", test = "t",dfs = "contain",struct = c("GEN","GEN"))

secondary.str.rom.rcs.model.ri=rma.mv(yi.rom,vi.rom,data = secondary.data.str,
                                    mods = ~ns(avg.rir,knots = secondary.str.k,Boundary.knots = secondary.str.b)+load.set+set.rep.equated+weeks+train.status,
                                    random = list(~1|study,~1|group,~1|obs),
                                    method = "REML", test = "t",dfs = "contain")

secondary.str.rom.rcs.model.rs=rma.mv(yi.rom,vi.rom,data = secondary.data.str,
                                    mods = ~ns(avg.rir,knots = secondary.str.k,Boundary.knots = secondary.str.b)+load.set+set.rep.equated+weeks+train.status,
                                    random = list(~avg.rir|study,~1|group,~1|obs),
                                    method = "REML", test = "t",dfs = "contain",struct = "GEN")

secondary.str.rom.rcs.model.rs2=rma.mv(yi.rom,vi.rom,data = secondary.data.str,
                                     mods = ~ns(avg.rir,knots = secondary.str.k,Boundary.knots = secondary.str.b)+load.set+set.rep.equated+weeks+train.status,
                                     random = list(~avg.rir|study,~avg.rir|group,~1|obs),
                                     method = "REML", test = "t",dfs = "contain",struct = c("GEN","GEN"))

#### compare g models
secondary.model.compare.str.g=bf_models(secondary.str.g.linear.model.ri,
                                      secondary.str.g.linear.model.rs,
                                      secondary.str.g.linear.model.rs2,
                                      secondary.str.g.log.model.ri,
                                      secondary.str.g.log.model.rs,
                                      secondary.str.g.log.model.rs2,
                                      secondary.str.g.quadratic.model.ri,
                                      secondary.str.g.quadratic.model.rs,
                                      secondary.str.g.quadratic.model.rs2,
                                      secondary.str.g.cubic.model.ri,
                                      secondary.str.g.cubic.model.rs,
                                      secondary.str.g.cubic.model.rs2,
                                      secondary.str.g.spline.model.ri,
                                      secondary.str.g.spline.model.rs,
                                      secondary.str.g.spline.model.rs2,
                                      secondary.str.g.rcs.model.ri,
                                      secondary.str.g.rcs.model.rs,
                                      secondary.str.g.rcs.model.rs2
)%>%as.matrix()%>%as.data.frame()%>%
  rename(linear.ri=1,
         linear.rs=2,
         linear.rs2=3,
         log.ri=4,
         log.rs=5,
         log.rs2=6,
         quadratic.ri=7,
         quadratic.rs=8,
         quadratic.rs2=9,
         cubic.ri=10,
         cubic.rs=11,
         cubic.rs2=12,
         spline.ri=13,
         spline.rs=14,
         spline.rs2=15,
         rcs.ri=16,
         rcs.rs=17,
         rcs.rs2=18)%>%
  mutate(model=c(
    "linear.ri",
    "linear.rs",
    "linear.rs2",
    "log.ri",
    "log.rs",
    "log.rs2",
    "quadratic.ri",
    "quadratic.rs",
    "quadratic.rs2",
    "cubic.ri",
    "cubic.rs",
    "cubic.rs2",
    "spline.ri",
    "spline.rs",
    "spline.rs2",
    "rcs.ri",
    "rcs.rs",
    "rcs.rs2"))%>%
  remove_rownames()%>%
  melt()%>%
  mutate(value=round(2*value,2))%>% 
  ggplot(aes(x=fct_inorder(variable),
             y=fct_relevel(model,
                           "rcs.rs2",
                           "rcs.rs",
                           "rcs.ri",
                           "spline.rs2",
                           "spline.rs",
                           "spline.ri",
                           "cubic.rs2",
                           "cubic.rs",
                           "cubic.ri",
                           "quadratic.rs2",
                           "quadratic.rs",
                           "quadratic.ri",
                           "log.rs2",
                           "log.rs",
                           "log.ri",
                           "linear.rs2",
                           "linear.rs",
                           "linear.ri")))+THEME+theme(legend.position = "right",
                                                      plot.subtitle = element_text(face = c("italic")))+
  labs(x="Numerator",
       y="Denominator",
       title = "Maximal Strength Model Comparision Matrix with <10 RIR (Standardized Mean Change)",
       subtitle = "2 x Log(Bayes Factor) Positive Values Favor Numerator",
       fill="2 x Log(BF)",
       caption = "Kass and Raferty (1995) Scale: -Inf to 0 = Negative; 0 to 2 = Weak; 2 to 6 = Positive; 6 to 10 = Strong; 10 to +Inf = Very Strong")+
  geom_tile(aes(fill=value))+
  geom_text(aes(label=value),color="black")+
  scale_fill_gradient2(low = "gray50",
                       high = "#c1272d",
                       mid="white",
                       midpoint = 0)+
  scale_x_discrete(labels=c("Linear\n(Intercept)",
                            "Linear\n(Slope)",
                            "Linear\n(Slopes)",
                            "Linear-log\n(Intercept)",
                            "Linear-log\n(Slope)",
                            "Linear-log\n(Slopes)",
                            "Quadratic\n(Intercept)",
                            "Quadratic\n(Slope)",
                            "Quadratic\n(Slopes)",
                            "Cubic\n(Intercept)",
                            "Cubic\n(Slope)",
                            "Cubic\n(Slopes)",
                            "Linear\nSpline\n(Intercept)",
                            "Linear\nSpline\n(Slope)",
                            "Linear\nSpline\n(Slopes)",
                            "Restricted\nCubic Spline\n(Intercept)",
                            "Restricted\nCubic Spline\n(Slope)",
                            "Restricted\nCubic Spline\n(Slopes)"))+
  scale_y_discrete(labels=rev(c("Linear\n(Intercept)",
                                "Linear\n(Slope)",
                                "Linear\n(Slopes)",
                                "Linear-log\n(Intercept)",
                                "Linear-log\n(Slope)",
                                "Linear-log\n(Slopes)",
                                "Quadratic\n(Intercept)",
                                "Quadratic\n(Slope)",
                                "Quadratic\n(Slopes)",
                                "Cubic\n(Intercept)",
                                "Cubic\n(Slope)",
                                "Cubic\n(Slopes)",
                                "Linear\nSpline\n(Intercept)",
                                "Linear\nSpline\n(Slope)",
                                "Linear\nSpline\n(Slopes)",
                                "Restricted\nCubic Spline\n(Intercept)",
                                "Restricted\nCubic Spline\n(Slope)", 
                                "Restricted\nCubic Spline\n(Slopes)"))) # log.ri best fit

ggsave(secondary.model.compare.str.g,device = "pdf",filename="secondary.model.compare.str.g.pdf",width = 14,height = 9,
       path = "/Users/zacrobinson/Documents/For Publication PTF Meta/Plots")

#### refit g model with robust variance estimates
secondary.str.g.best.model.prep=rma.mv(yi,vi,data = secondary.data.str,
                                     mods = ~ log1p(avg.rir)+load.set+set.rep.equated+weeks+train.status,
                                     random = list(~1|study,~1|group,~1|obs),
                                     method = "REML", test = "t",dfs = "contain")

secondary.str.g.best.model=robust.rma.mv(secondary.str.g.best.model.prep,
                                       cluster = secondary.data.str$study)

#### I2 and R2 for g model
secondary.str.g.best.model.r2=r2_ml(secondary.str.g.best.model)*100
save(secondary.str.g.best.model.r2,file = "secondary.str.g.best.model.r2.RData")
secondary.str.g.best.model.i2=i2_ml(secondary.str.g.best.model)
save(secondary.str.g.best.model.i2,file = "secondary.str.g.best.model.i2.RData")

#### g model output
secondary.str.g.best.output<-subset(as.data.frame(coef(summary(secondary.str.g.best.model))),select = c(-pval,-tval))
secondary.str.g.best.output.prep1<-qt(.975,secondary.str.g.best.output$df)
secondary.str.g.best.output.prep2<-sqrt((secondary.str.g.best.output$se^2)+sum(secondary.str.g.best.model$sigma2))
secondary.str.g.best.output$lower.PI<-secondary.str.g.best.output$estimate-(secondary.str.g.best.output.prep1*secondary.str.g.best.output.prep2)
secondary.str.g.best.output$upper.PI<-secondary.str.g.best.output$estimate+(secondary.str.g.best.output.prep1*secondary.str.g.best.output.prep2)
secondary.str.g.best.output$Parameter<-c("Intercept",
                                       "Log<sub>+1</sub>(Estimated RIR)",
                                       "Load",
                                       "Volume Equated [Rep]",
                                       "Volume Equated [Both]",
                                       "Weeks",
                                       "Training Status [Untrained]")

secondary.str.g.best.output.table<-relocate(secondary.str.g.best.output,Parameter)%>%
  gt()%>%
  cols_label(ci.lb="Lower Bound",
             ci.ub="Upper Bound",
             lower.PI="Lower Bound",
             upper.PI="Upper Bound",
             se="Standard Error",
             estimate="Estimate")%>%
  fmt_markdown(columns = 1)%>%
  tab_header(title = "Coefficients of Linear-log Multilevel Model for Maximal Strength with <10 RIR (Standardized Mean Change)")%>%
  tab_spanner(columns = 5:6,label = "95% Confidence Interval")%>%
  tab_spanner(columns = 7:8,label = "95% Predicition Interval")%>%
  opt_stylize(style = 1,color = "gray")%>%
  tab_options(column_labels.background.color ="gray97",column_labels.font.weight = "bold")%>%
  tab_style(style = cell_text(style = "italic"),locations = cells_row_groups())%>%
  tab_style(style = cell_text(align = "center"),locations = cells_column_labels())%>%
  tab_style(style = cell_text(weight  = "bold"),locations = cells_title())%>%
  tab_style(style = cell_text(align = "center"),locations = cells_body(columns=-1))%>%
  cols_width(-df~px(133))%>%
  cols_width(Parameter~px(200))

gtsave(secondary.str.g.best.output.table,filename = "secondary.str.g.best.output.table.png",
       path = "/Users/zacrobinson/Documents/For Publication PTF Meta/Models")

#### compare rom models
secondary.model.compare.str.rom=bf_models(secondary.str.rom.linear.model.ri,
                                        secondary.str.rom.linear.model.rs,
                                        secondary.str.rom.linear.model.rs2,
                                        secondary.str.rom.log.model.ri,
                                        secondary.str.rom.log.model.rs,
                                        secondary.str.rom.log.model.rs2,
                                        secondary.str.rom.quadratic.model.ri,
                                        secondary.str.rom.quadratic.model.rs,
                                        secondary.str.rom.quadratic.model.rs2,
                                        secondary.str.rom.cubic.model.ri,
                                        secondary.str.rom.cubic.model.rs,
                                        secondary.str.rom.cubic.model.rs2,
                                        secondary.str.rom.spline.model.ri,
                                        secondary.str.rom.spline.model.rs,
                                        secondary.str.rom.spline.model.rs2,
                                        secondary.str.rom.rcs.model.ri,
                                        secondary.str.rom.rcs.model.rs,
                                        secondary.str.rom.rcs.model.rs2
)%>%as.matrix()%>%as.data.frame()%>%
  rename(linear.ri=1,
         linear.rs=2,
         linear.rs2=3,
         log.ri=4,
         log.rs=5,
         log.rs2=6,
         quadratic.ri=7,
         quadratic.rs=8,
         quadratic.rs2=9,
         cubic.ri=10,
         cubic.rs=11,
         cubic.rs2=12,
         spline.ri=13,
         spline.rs=14,
         spline.rs2=15,
         rcs.ri=16,
         rcs.rs=17,
         rcs.rs2=18)%>%
  mutate(model=c(
    "linear.ri",
    "linear.rs",
    "linear.rs2",
    "log.ri",
    "log.rs",
    "log.rs2",
    "quadratic.ri",
    "quadratic.rs",
    "quadratic.rs2",
    "cubic.ri",
    "cubic.rs",
    "cubic.rs2",
    "spline.ri",
    "spline.rs",
    "spline.rs2",
    "rcs.ri",
    "rcs.rs",
    "rcs.rs2"))%>%
  remove_rownames()%>%
  melt()%>%
  mutate(value=round(2*value,2))%>% 
  ggplot(aes(x=fct_inorder(variable),
             y=fct_relevel(model,
                           "rcs.rs2",
                           "rcs.rs",
                           "rcs.ri",
                           "spline.rs2",
                           "spline.rs",
                           "spline.ri",
                           "cubic.rs2",
                           "cubic.rs",
                           "cubic.ri",
                           "quadratic.rs2",
                           "quadratic.rs",
                           "quadratic.ri",
                           "log.rs2",
                           "log.rs",
                           "log.ri",
                           "linear.rs2",
                           "linear.rs",
                           "linear.ri")))+THEME+theme(legend.position = "right",
                                                      plot.subtitle = element_text(face = c("italic")))+
  labs(x="Numerator",
       y="Denominator",
       title = "Maximal Strength Model Comparision Matrix with <10 RIR (Response Ratio)",
       subtitle = "2 x Log(Bayes Factor) Positive Values Favor Numerator",
       fill="2 x Log(BF)",
       caption = "Kass and Raferty (1995) Scale: -Inf to 0 = Negative; 0 to 2 = Weak; 2 to 6 = Positive; 6 to 10 = Strong; 10 to +Inf = Very Strong")+
  geom_tile(aes(fill=value))+
  geom_text(aes(label=value),color="black")+
  scale_fill_gradient2(low = "gray50",
                       high = "#c1272d",
                       mid="white",
                       midpoint = 0)+
  scale_x_discrete(labels=c("Linear\n(Intercept)",
                            "Linear\n(Slope)",
                            "Linear\n(Slopes)",
                            "Linear-log\n(Intercept)",
                            "Linear-log\n(Slope)",
                            "Linear-log\n(Slopes)",
                            "Quadratic\n(Intercept)",
                            "Quadratic\n(Slope)",
                            "Quadratic\n(Slopes)",
                            "Cubic\n(Intercept)",
                            "Cubic\n(Slope)",
                            "Cubic\n(Slopes)",
                            "Linear\nSpline\n(Intercept)",
                            "Linear\nSpline\n(Slope)",
                            "Linear\nSpline\n(Slopes)",
                            "Restricted\nCubic Spline\n(Intercept)",
                            "Restricted\nCubic Spline\n(Slope)",
                            "Restricted\nCubic Spline\n(Slopes)"))+
  scale_y_discrete(labels=rev(c("Linear\n(Intercept)",
                                "Linear\n(Slope)",
                                "Linear\n(Slopes)",
                                "Linear-log\n(Intercept)",
                                "Linear-log\n(Slope)",
                                "Linear-log\n(Slopes)",
                                "Quadratic\n(Intercept)",
                                "Quadratic\n(Slope)",
                                "Quadratic\n(Slopes)",
                                "Cubic\n(Intercept)",
                                "Cubic\n(Slope)",
                                "Cubic\n(Slopes)",
                                "Linear\nSpline\n(Intercept)",
                                "Linear\nSpline\n(Slope)",
                                "Linear\nSpline\n(Slopes)",
                                "Restricted\nCubic Spline\n(Intercept)",
                                "Restricted\nCubic Spline\n(Slope)", 
                                "Restricted\nCubic Spline\n(Slopes)"))) # linear.ri best fit

ggsave(secondary.model.compare.str.rom,device = "pdf",filename="secondary.model.compare.str.rom.pdf",width = 14,height = 9,
       path = "/Users/zacrobinson/Documents/For Publication PTF Meta/Plots")

#### refit rom model with robust variance estimates
secondary.str.rom.best.model.prep=rma.mv(yi.rom,vi.rom,data = secondary.data.str,
                                       mods = ~avg.rir+load.set+set.rep.equated+weeks+train.status,
                                       random = list(~1|study,~1|group,~1|obs),
                                       method = "REML", test = "t",dfs = "contain")

secondary.str.rom.best.model=robust.rma.mv(secondary.str.rom.best.model.prep,
                                         cluster = secondary.data.str$study)

#### I2 and R2 for rom model
secondary.str.rom.best.model.r2=r2_ml(secondary.str.rom.best.model)*100
save(secondary.str.rom.best.model.r2,file = "secondary.str.rom.best.model.r2.RData")
secondary.str.rom.best.model.i2=i2_ml(secondary.str.rom.best.model)
save(secondary.str.rom.best.model.i2,file = "secondary.str.rom.best.model.i2.RData")

#### rom model output
secondary.str.rom.best.output<-subset(as.data.frame(coef(summary(secondary.str.rom.best.model))),select = c(-pval,-tval))
secondary.str.rom.best.output.prep1<-qt(.975,secondary.str.rom.best.output$df)
secondary.str.rom.best.output.prep2<-sqrt((secondary.str.rom.best.output$se^2)+sum(secondary.str.rom.best.model$sigma2))
secondary.str.rom.best.output$lower.PI<-secondary.str.rom.best.output$estimate-(secondary.str.rom.best.output.prep1*secondary.str.rom.best.output.prep2)
secondary.str.rom.best.output$upper.PI<-secondary.str.rom.best.output$estimate+(secondary.str.rom.best.output.prep1*secondary.str.rom.best.output.prep2)
secondary.str.rom.best.output$Parameter<-c("Intercept",
                                         "Estimated RIR",
                                         "Load",
                                         "Volume Equated [Rep]",
                                         "Volume Equated [Both]",
                                         "Weeks",
                                         "Training Status [Untrained]")

secondary.str.rom.best.output.table<-relocate(secondary.str.rom.best.output,Parameter)%>%
  gt()%>%
  cols_label(ci.lb="Lower Bound",
             ci.ub="Upper Bound",
             lower.PI="Lower Bound",
             upper.PI="Upper Bound",
             se="Standard Error",
             estimate="Estimate")%>%
  fmt_markdown(columns = 1)%>%
  tab_header(title = "Coefficients of Linear Multilevel Model for Maximal Strength with <10 RIR (Response Ratio)")%>%
  tab_spanner(columns = 5:6,label = "95% Confidence Interval")%>%
  tab_spanner(columns = 7:8,label = "95% Predicition Interval")%>%
  opt_stylize(style = 1,color = "gray")%>%
  tab_options(column_labels.background.color ="gray97",column_labels.font.weight = "bold")%>%
  tab_style(style = cell_text(style = "italic"),locations = cells_row_groups())%>%
  tab_style(style = cell_text(align = "center"),locations = cells_column_labels())%>%
  tab_style(style = cell_text(weight  = "bold"),locations = cells_title())%>%
  tab_style(style = cell_text(align = "center"),locations = cells_body(columns=-1))%>%
  cols_width(-df~px(133))%>%
  cols_width(Parameter~px(200))

gtsave(secondary.str.rom.best.output.table,filename = "secondary.str.rom.best.output.table.png",
       path = "/Users/zacrobinson/Documents/For Publication PTF Meta/Models")

### hypertrophy
#### linear
secondary.hyp.g.linear.model.ri=rma.mv(yi,vi,data = secondary.data.hyp,
                                     mods = ~ avg.rir+load.set+set.rep.equated+weeks+train.status,
                                     random = list(~1|study,~1|group,~1|obs),
                                     method = "REML", test = "t",dfs = "contain")

secondary.hyp.g.linear.model.rs=rma.mv(yi,vi,data = secondary.data.hyp,
                                     mods = ~ avg.rir+load.set+set.rep.equated+weeks+train.status,
                                     random = list(~avg.rir|study,~1|group,~1|obs),
                                     method = "REML", test = "t",dfs = "contain",struct = "GEN")

secondary.hyp.g.linear.model.rs2=rma.mv(yi,vi,data = secondary.data.hyp,
                                      mods = ~ avg.rir+load.set+set.rep.equated+weeks+train.status,
                                      random = list(~avg.rir|study,~avg.rir|group,~1|obs),
                                      method = "REML", test = "t",dfs = "contain",struct = c("GEN","GEN"))

secondary.hyp.rom.linear.model.ri=rma.mv(yi.rom,vi.rom,data = secondary.data.hyp,
                                       mods = ~ avg.rir+load.set+set.rep.equated+weeks+train.status,
                                       random = list(~1|study,~1|group,~1|obs),
                                       method = "REML", test = "t",dfs = "contain")

secondary.hyp.rom.linear.model.rs=rma.mv(yi.rom,vi.rom,data = secondary.data.hyp,
                                       mods = ~ avg.rir+load.set+set.rep.equated+weeks+train.status,
                                       random = list(~avg.rir|study,~1|group,~1|obs),
                                       method = "REML", test = "t",dfs = "contain",struct = "GEN")

secondary.hyp.rom.linear.model.rs2=rma.mv(yi.rom,vi.rom,data = secondary.data.hyp,
                                        mods = ~ avg.rir+load.set+set.rep.equated+weeks+train.status,
                                        random = list(~avg.rir|study,~avg.rir|group,~1|obs),
                                        method = "REML", test = "t",dfs = "contain",struct = c("GEN","GEN"))

#### linear-log
secondary.hyp.g.log.model.ri=rma.mv(yi,vi,data = secondary.data.hyp,
                                  mods = ~ log1p(avg.rir)+load.set+set.rep.equated+weeks+train.status,
                                  random = list(~1|study,~1|group,~1|obs),
                                  method = "REML", test = "t",dfs = "contain")

secondary.hyp.g.log.model.rs=rma.mv(yi,vi,data = secondary.data.hyp,
                                  mods = ~ log1p(avg.rir)+load.set+set.rep.equated+weeks+train.status,
                                  random = list(~avg.rir|study,~1|group,~1|obs),
                                  method = "REML", test = "t",dfs = "contain",struct = "GEN")

secondary.hyp.g.log.model.rs2=rma.mv(yi,vi,data = secondary.data.hyp,
                                   mods = ~ log1p(avg.rir)+load.set+set.rep.equated+weeks+train.status,
                                   random = list(~avg.rir|study,~avg.rir|group,~1|obs),
                                   method = "REML", test = "t",dfs = "contain",struct = c("GEN","GEN"))

secondary.hyp.rom.log.model.ri=rma.mv(yi.rom,vi.rom,data = secondary.data.hyp,
                                    mods = ~ log1p(avg.rir)+load.set+set.rep.equated+weeks+train.status,
                                    random = list(~1|study,~1|group,~1|obs),
                                    method = "REML", test = "t",dfs = "contain")

secondary.hyp.rom.log.model.rs=rma.mv(yi.rom,vi.rom,data = secondary.data.hyp,
                                    mods = ~ log1p(avg.rir)+load.set+set.rep.equated+weeks+train.status,
                                    random = list(~avg.rir|study,~1|group,~1|obs),
                                    method = "REML", test = "t",dfs = "contain",struct = "GEN")

secondary.hyp.rom.log.model.rs2=rma.mv(yi.rom,vi.rom,data = secondary.data.hyp,
                                     mods = ~ log1p(avg.rir)+load.set+set.rep.equated+weeks+train.status,
                                     random = list(~avg.rir|study,~avg.rir|group,~1|obs),
                                     method = "REML", test = "t",dfs = "contain",struct = c("GEN","GEN"))

#### quadratic
secondary.hyp.g.quadratic.model.ri=rma.mv(yi,vi,data = secondary.data.hyp,
                                        mods = ~ poly(avg.rir,2,raw = TRUE)+load.set+set.rep.equated+weeks+train.status,
                                        random = list(~1|study,~1|group,~1|obs),
                                        method = "REML", test = "t",dfs = "contain")

secondary.hyp.g.quadratic.model.rs=rma.mv(yi,vi,data = secondary.data.hyp,
                                        mods = ~ poly(avg.rir,2,raw = TRUE)+load.set+set.rep.equated+weeks+train.status,
                                        random = list(~avg.rir|study,~1|group,~1|obs),
                                        method = "REML", test = "t",dfs = "contain",struct = "GEN")

secondary.hyp.g.quadratic.model.rs2=rma.mv(yi,vi,data = secondary.data.hyp,
                                         mods = ~ poly(avg.rir,2,raw = TRUE)+load.set+set.rep.equated+weeks+train.status,
                                         random = list(~avg.rir|study,~avg.rir|group,~1|obs),
                                         method = "REML", test = "t",dfs = "contain",struct = c("GEN","GEN"))

secondary.hyp.rom.quadratic.model.ri=rma.mv(yi.rom,vi.rom,data = secondary.data.hyp,
                                          mods = ~ poly(avg.rir,2,raw = TRUE)+load.set+set.rep.equated+weeks+train.status,
                                          random = list(~1|study,~1|group,~1|obs),
                                          method = "REML", test = "t",dfs = "contain")

secondary.hyp.rom.quadratic.model.rs=rma.mv(yi.rom,vi.rom,data = secondary.data.hyp,
                                          mods = ~ poly(avg.rir,2,raw = TRUE)+load.set+set.rep.equated+weeks+train.status,
                                          random = list(~avg.rir|study,~1|group,~1|obs),
                                          method = "REML", test = "t",dfs = "contain",struct = "GEN")

secondary.hyp.rom.quadratic.model.rs2=rma.mv(yi.rom,vi.rom,data = secondary.data.hyp,
                                           mods = ~ poly(avg.rir,2,raw = TRUE)+load.set+set.rep.equated+weeks+train.status,
                                           random = list(~avg.rir|study,~avg.rir|group,~1|obs),
                                           method = "REML", test = "t",dfs = "contain",struct = c("GEN","GEN"))

#### cubic
secondary.hyp.g.cubic.model.ri=rma.mv(yi,vi,data = secondary.data.hyp,
                                    mods = ~ poly(avg.rir,3,raw = TRUE)+load.set+set.rep.equated+weeks+train.status,
                                    random = list(~1|study,~1|group,~1|obs),
                                    method = "REML", test = "t",dfs = "contain")

secondary.hyp.g.cubic.model.rs=rma.mv(yi,vi,data = secondary.data.hyp,
                                    mods = ~ poly(avg.rir,3,raw = TRUE)+load.set+set.rep.equated+weeks+train.status,
                                    random = list(~avg.rir|study,~1|group,~1|obs),
                                    method = "REML", test = "t",dfs = "contain",struct = "GEN")

secondary.hyp.g.cubic.model.rs2=rma.mv(yi,vi,data = secondary.data.hyp,
                                     mods = ~ poly(avg.rir,3,raw = TRUE)+load.set+set.rep.equated+weeks+train.status,
                                     random = list(~avg.rir|study,~avg.rir|group,~1|obs),
                                     method = "REML", test = "t",dfs = "contain",struct = c("GEN","GEN"))

secondary.hyp.rom.cubic.model.ri=rma.mv(yi.rom,vi.rom,data = secondary.data.hyp,
                                      mods = ~ poly(avg.rir,3,raw = TRUE)+load.set+set.rep.equated+weeks+train.status,
                                      random = list(~1|study,~1|group,~1|obs),
                                      method = "REML", test = "t",dfs = "contain")

secondary.hyp.rom.cubic.model.rs=rma.mv(yi.rom,vi.rom,data = secondary.data.hyp,
                                      mods = ~ poly(avg.rir,3,raw = TRUE)+load.set+set.rep.equated+weeks+train.status,
                                      random = list(~avg.rir|study,~1|group,~1|obs),
                                      method = "REML", test = "t",dfs = "contain",struct = "GEN")

secondary.hyp.rom.cubic.model.rs2=rma.mv(yi.rom,vi.rom,data = secondary.data.hyp,
                                       mods = ~ poly(avg.rir,3,raw = TRUE)+load.set+set.rep.equated+weeks+train.status,
                                       random = list(~avg.rir|study,~avg.rir|group,~1|obs),
                                       method = "REML", test = "t",dfs = "contain",struct = c("GEN","GEN"))

#### linear-spline
secondary.hyp.g.spline.model.ri=rma.mv(yi,vi,data = secondary.data.hyp,
                                     mods = ~ bs(spline.rir,degree = 1,knots = 0)+load.set+set.rep.equated+weeks+train.status,
                                     random = list(~1|study,~1|group,~1|obs),
                                     method = "REML", test = "t",dfs = "contain")

secondary.hyp.g.spline.model.rs=rma.mv(yi,vi,data = secondary.data.hyp,
                                     mods = ~ bs(spline.rir,degree = 1,knots = 0)+load.set+set.rep.equated+weeks+train.status,
                                     random = list(~spline.rir|study,~1|group,~1|obs),
                                     method = "REML", test = "t",dfs = "contain",struct = "GEN")

secondary.hyp.g.spline.model.rs2=rma.mv(yi,vi,data = secondary.data.hyp,
                                      mods = ~ bs(spline.rir,degree = 1,knots = 0)+load.set+set.rep.equated+weeks+train.status,
                                      random = list(~spline.rir|study,~spline.rir|group,~1|obs),
                                      method = "REML", test = "t",dfs = "contain",struct = c("GEN","GEN"))

secondary.hyp.rom.spline.model.ri=rma.mv(yi.rom,vi.rom,data = secondary.data.hyp,
                                       mods = ~ bs(spline.rir,degree = 1,knots = 0)+load.set+set.rep.equated+weeks+train.status,
                                       random = list(~1|study,~1|group,~1|obs),
                                       method = "REML", test = "t",dfs = "contain")

secondary.hyp.rom.spline.model.rs=rma.mv(yi.rom,vi.rom,data = secondary.data.hyp,
                                       mods = ~ bs(spline.rir,degree = 1,knots = 0)+load.set+set.rep.equated+weeks+train.status,
                                       random = list(~spline.rir|study,~1|group,~1|obs),
                                       method = "REML", test = "t",dfs = "contain",struct = "GEN")

secondary.hyp.rom.spline.model.rs2=rma.mv(yi.rom,vi.rom,data = secondary.data.hyp,
                                        mods = ~ bs(spline.rir,degree = 1,knots = 0)+load.set+set.rep.equated+weeks+train.status,
                                        random = list(~spline.rir|study,~spline.rir|group,~1|obs),
                                        method = "REML", test = "t",dfs = "contain",struct = c("GEN","GEN"))

#### cubic-spline
secondary.hyp.g.rcs.model.ri=rma.mv(yi,vi,data = secondary.data.hyp,
                                  mods = ~ ns(avg.rir,knots = secondary.hyp.k,Boundary.knots = secondary.hyp.b)+load.set+set.rep.equated+weeks+train.status,
                                  random = list(~1|study,~1|group,~1|obs),
                                  method = "REML", test = "t",dfs = "contain")

secondary.hyp.g.rcs.model.rs=rma.mv(yi,vi,data = secondary.data.hyp,
                                  mods = ~ ns(avg.rir,knots = secondary.hyp.k,Boundary.knots = secondary.hyp.b)+load.set+set.rep.equated+weeks+train.status,
                                  random = list(~avg.rir|study,~1|group,~1|obs),
                                  method = "REML", test = "t",dfs = "contain",struct = "GEN")

secondary.hyp.g.rcs.model.rs2=rma.mv(yi,vi,data = secondary.data.hyp,
                                   mods = ~ ns(avg.rir,knots = secondary.hyp.k,Boundary.knots = secondary.hyp.b)+load.set+set.rep.equated+weeks+train.status,
                                   random = list(~avg.rir|study,~avg.rir|group,~1|obs),
                                   method = "REML", test = "t",dfs = "contain",struct = c("GEN","GEN"))

secondary.hyp.rom.rcs.model.ri=rma.mv(yi.rom,vi.rom,data = secondary.data.hyp,
                                    mods = ~ns(avg.rir,knots = secondary.hyp.k,Boundary.knots = secondary.hyp.b)+load.set+set.rep.equated+weeks+train.status,
                                    random = list(~1|study,~1|group,~1|obs),
                                    method = "REML", test = "t",dfs = "contain")

secondary.hyp.rom.rcs.model.rs=rma.mv(yi.rom,vi.rom,data = secondary.data.hyp,
                                    mods = ~ns(avg.rir,knots = secondary.hyp.k,Boundary.knots = secondary.hyp.b)+load.set+set.rep.equated+weeks+train.status,
                                    random = list(~avg.rir|study,~1|group,~1|obs),
                                    method = "REML", test = "t",dfs = "contain",struct = "GEN",)

secondary.hyp.rom.rcs.model.rs2=rma.mv(yi.rom,vi.rom,data = secondary.data.hyp,
                                     mods = ~ns(avg.rir,knots = secondary.hyp.k,Boundary.knots = secondary.hyp.b)+load.set+set.rep.equated+weeks+train.status,
                                     random = list(~avg.rir|study,~avg.rir|group,~1|obs),
                                     method = "REML", test = "t",dfs = "contain",struct = c("GEN","GEN"))

#### compare models
secondary.model.compare.hyp.g=bf_models(secondary.hyp.g.linear.model.ri,
                                      secondary.hyp.g.linear.model.rs,
                                      secondary.hyp.g.linear.model.rs2,
                                      secondary.hyp.g.log.model.ri,
                                      secondary.hyp.g.log.model.rs,
                                      secondary.hyp.g.log.model.rs2,
                                      secondary.hyp.g.quadratic.model.ri,
                                      secondary.hyp.g.quadratic.model.rs,
                                      secondary.hyp.g.quadratic.model.rs2,
                                      secondary.hyp.g.cubic.model.ri,
                                      secondary.hyp.g.cubic.model.rs,
                                      secondary.hyp.g.cubic.model.rs2,
                                      secondary.hyp.g.spline.model.ri,
                                      secondary.hyp.g.spline.model.rs,
                                      secondary.hyp.g.spline.model.rs2,
                                      secondary.hyp.g.rcs.model.ri,
                                      secondary.hyp.g.rcs.model.rs,
                                      secondary.hyp.g.rcs.model.rs2
)%>%as.matrix()%>%as.data.frame()%>%
  rename(linear.ri=1,
         linear.rs=2,
         linear.rs2=3,
         log.ri=4,
         log.rs=5,
         log.rs2=6,
         quadratic.ri=7,
         quadratic.rs=8,
         quadratic.rs2=9,
         cubic.ri=10,
         cubic.rs=11,
         cubic.rs2=12,
         spline.ri=13,
         spline.rs=14,
         spline.rs2=15,
         rcs.ri=16,
         rcs.rs=17,
         rcs.rs2=18)%>%
  mutate(model=c(
    "linear.ri",
    "linear.rs",
    "linear.rs2",
    "log.ri",
    "log.rs",
    "log.rs2",
    "quadratic.ri",
    "quadratic.rs",
    "quadratic.rs2",
    "cubic.ri",
    "cubic.rs",
    "cubic.rs2",
    "spline.ri",
    "spline.rs",
    "spline.rs2",
    "rcs.ri",
    "rcs.rs",
    "rcs.rs2"))%>%
  remove_rownames()%>%
  melt()%>%
  mutate(value=round(2*value,2))%>% 
  ggplot(aes(x=fct_inorder(variable),
             y=fct_relevel(model,
                           "rcs.rs2",
                           "rcs.rs",
                           "rcs.ri",
                           "spline.rs2",
                           "spline.rs",
                           "spline.ri",
                           "cubic.rs2",
                           "cubic.rs",
                           "cubic.ri",
                           "quadratic.rs2",
                           "quadratic.rs",
                           "quadratic.ri",
                           "log.rs2",
                           "log.rs",
                           "log.ri",
                           "linear.rs2",
                           "linear.rs",
                           "linear.ri")))+THEME+theme(legend.position = "right",
                                                      plot.subtitle = element_text(face = c("italic")))+
  labs(x="Numerator",
       y="Denominator",
       title = "Muscle Hypertrophy Model Comparision Matrix with <10 RIR (Standardized Mean Change)",
       subtitle = "2 x Log(Bayes Factor) Positive Values Favor Numerator",
       fill="2 x Log(BF)",
       caption = "Kass and Raferty (1995) Scale: -Inf to 0 = Negative; 0 to 2 = Weak; 2 to 6 = Positive; 6 to 10 = Strong; 10 to +Inf = Very Strong")+
  geom_tile(aes(fill=value))+
  geom_text(aes(label=value),color="black")+
  scale_fill_gradient2(low = "gray50",
                       high = "#0000a7",
                       mid="white",
                       midpoint = 0)+
  scale_x_discrete(labels=c("Linear\n(Intercept)",
                            "Linear\n(Slope)",
                            "Linear\n(Slopes)",
                            "Linear-log\n(Intercept)",
                            "Linear-log\n(Slope)",
                            "Linear-log\n(Slopes)",
                            "Quadratic\n(Intercept)",
                            "Quadratic\n(Slope)",
                            "Quadratic\n(Slopes)",
                            "Cubic\n(Intercept)",
                            "Cubic\n(Slope)",
                            "Cubic\n(Slopes)",
                            "Linear\nSpline\n(Intercept)",
                            "Linear\nSpline\n(Slope)",
                            "Linear\nSpline\n(Slopes)",
                            "Restricted\nCubic Spline\n(Intercept)",
                            "Restricted\nCubic Spline\n(Slope)",
                            "Restricted\nCubic Spline\n(Slopes)"))+
  scale_y_discrete(labels=rev(c("Linear\n(Intercept)",
                                "Linear\n(Slope)",
                                "Linear\n(Slopes)",
                                "Linear-log\n(Intercept)",
                                "Linear-log\n(Slope)",
                                "Linear-log\n(Slopes)",
                                "Quadratic\n(Intercept)",
                                "Quadratic\n(Slope)",
                                "Quadratic\n(Slopes)",
                                "Cubic\n(Intercept)",
                                "Cubic\n(Slope)",
                                "Cubic\n(Slopes)",
                                "Linear\nSpline\n(Intercept)",
                                "Linear\nSpline\n(Slope)",
                                "Linear\nSpline\n(Slopes)",
                                "Restricted\nCubic Spline\n(Intercept)",
                                "Restricted\nCubic Spline\n(Slope)", 
                                "Restricted\nCubic Spline\n(Slopes)"))) # linear.ri best fit

ggsave(secondary.model.compare.hyp.g,device = "pdf",filename="secondary.model.compare.hyp.g.pdf",width = 14,height = 9,
       path = "/Users/zacrobinson/Documents/For Publication PTF Meta/Plots")

#### refit g model with robust variance estimates
secondary.hyp.g.best.model.prep=rma.mv(yi,vi,data = secondary.data.hyp,
                                     mods = ~ avg.rir+load.set+set.rep.equated+weeks+train.status,
                                     random = list(~1|study,~1|group,~1|obs),
                                     method = "REML", test = "t",dfs = "contain")

secondary.hyp.g.best.model=robust.rma.mv(secondary.hyp.g.best.model.prep,
                                       cluster = secondary.data.hyp$study)
#### I2 and R2 for g model
secondary.hyp.g.best.model.r2=r2_ml(secondary.hyp.g.best.model)*100
save(secondary.hyp.g.best.model.r2,file = "secondary.hyp.g.best.model.r2.RData")
secondary.hyp.g.best.model.i2=i2_ml(secondary.hyp.g.best.model)
save(secondary.hyp.g.best.model.i2,file = "secondary.hyp.g.best.model.i2.RData")

#### g model output
secondary.hyp.g.best.output<-subset(as.data.frame(coef(summary(secondary.hyp.g.best.model))),select = c(-pval,-tval))
secondary.hyp.g.best.output.prep1<-qt(.975,secondary.hyp.g.best.output$df)
secondary.hyp.g.best.output.prep2<-sqrt((secondary.hyp.g.best.output$se^2)+sum(secondary.hyp.g.best.model$sigma2))
secondary.hyp.g.best.output$lower.PI<-secondary.hyp.g.best.output$estimate-(secondary.hyp.g.best.output.prep1*secondary.hyp.g.best.output.prep2)
secondary.hyp.g.best.output$upper.PI<-secondary.hyp.g.best.output$estimate+(secondary.hyp.g.best.output.prep1*secondary.hyp.g.best.output.prep2)
secondary.hyp.g.best.output$Parameter<-c("Intercept",
                                       "Estimated RIR",
                                       "Load",
                                       "Volume Equated [Rep]",
                                       "Volume Equated [Both]",
                                       "Weeks",
                                       "Training Status [Untrained]")

secondary.hyp.g.best.output.table<-relocate(secondary.hyp.g.best.output,Parameter)%>%
  gt()%>%
  cols_label(ci.lb="Lower Bound",
             ci.ub="Upper Bound",
             lower.PI="Lower Bound",
             upper.PI="Upper Bound",
             se="Standard Error",
             estimate="Estimate")%>%
  fmt_markdown(columns = 1)%>%
  tab_header(title = "Coefficients of Linear Multilevel Model for Muscle Hypertrophy with <10 RIR (Standardized Mean Change)")%>%
  tab_spanner(columns = 5:6,label = "95% Confidence Interval")%>%
  tab_spanner(columns = 7:8,label = "95% Predicition Interval")%>%
  opt_stylize(style = 1,color = "gray")%>%
  tab_options(column_labels.background.color ="gray97",column_labels.font.weight = "bold")%>%
  tab_style(style = cell_text(style = "italic"),locations = cells_row_groups())%>%
  tab_style(style = cell_text(align = "center"),locations = cells_column_labels())%>%
  tab_style(style = cell_text(weight  = "bold"),locations = cells_title())%>%
  tab_style(style = cell_text(align = "center"),locations = cells_body(columns=-1))%>%
  cols_width(-df~px(133))%>%
  cols_width(Parameter~px(200))

gtsave(secondary.hyp.g.best.output.table,filename = "secondary.hyp.g.best.output.table.png",
       path = "/Users/zacrobinson/Documents/For Publication PTF Meta/Models")

#### compare rom models
secondary.model.compare.hyp.rom=bf_models(secondary.hyp.rom.linear.model.ri,
                                        secondary.hyp.rom.linear.model.rs,
                                        secondary.hyp.rom.linear.model.rs2,
                                        secondary.hyp.rom.log.model.ri,
                                        secondary.hyp.rom.log.model.rs,
                                        secondary.hyp.rom.log.model.rs2,
                                        secondary.hyp.rom.quadratic.model.ri,
                                        secondary.hyp.rom.quadratic.model.rs,
                                        secondary.hyp.rom.quadratic.model.rs2,
                                        secondary.hyp.rom.cubic.model.ri,
                                        secondary.hyp.rom.cubic.model.rs,
                                        secondary.hyp.rom.cubic.model.rs2,
                                        secondary.hyp.rom.spline.model.ri,
                                        secondary.hyp.rom.spline.model.rs,
                                        secondary.hyp.rom.spline.model.rs2,
                                        secondary.hyp.rom.rcs.model.ri,
                                        secondary.hyp.rom.rcs.model.rs,
                                        secondary.hyp.rom.rcs.model.rs2
)%>%as.matrix()%>%as.data.frame()%>%
  rename(linear.ri=1,
         linear.rs=2,
         linear.rs2=3,
         log.ri=4,
         log.rs=5,
         log.rs2=6,
         quadratic.ri=7,
         quadratic.rs=8,
         quadratic.rs2=9,
         cubic.ri=10,
         cubic.rs=11,
         cubic.rs2=12,
         spline.ri=13,
         spline.rs=14,
         spline.rs2=15,
         rcs.ri=16,
         rcs.rs=17,
         rcs.rs2=18)%>%
  mutate(model=c(
    "linear.ri",
    "linear.rs",
    "linear.rs2",
    "log.ri",
    "log.rs",
    "log.rs2",
    "quadratic.ri",
    "quadratic.rs",
    "quadratic.rs2",
    "cubic.ri",
    "cubic.rs",
    "cubic.rs2",
    "spline.ri",
    "spline.rs",
    "spline.rs2",
    "rcs.ri",
    "rcs.rs",
    "rcs.rs2"))%>%
  remove_rownames()%>%
  melt()%>%
  mutate(value=round(2*value,2))%>% 
  ggplot(aes(x=fct_inorder(variable),
             y=fct_relevel(model,
                           "rcs.rs2",
                           "rcs.rs",
                           "rcs.ri",
                           "spline.rs2",
                           "spline.rs",
                           "spline.ri",
                           "cubic.rs2",
                           "cubic.rs",
                           "cubic.ri",
                           "quadratic.rs2",
                           "quadratic.rs",
                           "quadratic.ri",
                           "log.rs2",
                           "log.rs",
                           "log.ri",
                           "linear.rs2",
                           "linear.rs",
                           "linear.ri")))+THEME+theme(legend.position = "right",
                                                      plot.subtitle = element_text(face = c("italic")))+
  labs(x="Numerator",
       y="Denominator",
       title = "Muscle Hypertrophy Model Comparision Matrix with <10 RIR (Response Ratio)",
       subtitle = "2 x Log(Bayes Factor) Positive Values Favor Numerator",
       fill="2 x Log(BF)",
       caption = "Kass and Raferty (1995) Scale: -Inf to 0 = Negative; 0 to 2 = Weak; 2 to 6 = Positive; 6 to 10 = Strong; 10 to +Inf = Very Strong")+
  geom_tile(aes(fill=value))+
  geom_text(aes(label=value),color="black")+
  scale_fill_gradient2(low = "gray50",
                       high = "#0000a7",
                       mid="white",
                       midpoint = 0)+
  scale_x_discrete(labels=c("Linear\n(Intercept)",
                            "Linear\n(Slope)",
                            "Linear\n(Slopes)",
                            "Linear-log\n(Intercept)",
                            "Linear-log\n(Slope)",
                            "Linear-log\n(Slopes)",
                            "Quadratic\n(Intercept)",
                            "Quadratic\n(Slope)",
                            "Quadratic\n(Slopes)",
                            "Cubic\n(Intercept)",
                            "Cubic\n(Slope)",
                            "Cubic\n(Slopes)",
                            "Linear\nSpline\n(Intercept)",
                            "Linear\nSpline\n(Slope)",
                            "Linear\nSpline\n(Slopes)",
                            "Restricted\nCubic Spline\n(Intercept)",
                            "Restricted\nCubic Spline\n(Slope)",
                            "Restricted\nCubic Spline\n(Slopes)"))+
  scale_y_discrete(labels=rev(c("Linear\n(Intercept)",
                                "Linear\n(Slope)",
                                "Linear\n(Slopes)",
                                "Linear-log\n(Intercept)",
                                "Linear-log\n(Slope)",
                                "Linear-log\n(Slopes)",
                                "Quadratic\n(Intercept)",
                                "Quadratic\n(Slope)",
                                "Quadratic\n(Slopes)",
                                "Cubic\n(Intercept)",
                                "Cubic\n(Slope)",
                                "Cubic\n(Slopes)",
                                "Linear\nSpline\n(Intercept)",
                                "Linear\nSpline\n(Slope)",
                                "Linear\nSpline\n(Slopes)",
                                "Restricted\nCubic Spline\n(Intercept)",
                                "Restricted\nCubic Spline\n(Slope)", 
                                "Restricted\nCubic Spline\n(Slopes)"))) # linear.ri best fit

ggsave(secondary.model.compare.hyp.rom,device = "pdf",filename="secondary.model.compare.hyp.rom.pdf",width = 14,height = 9,
       path = "/Users/zacrobinson/Documents/For Publication PTF Meta/Plots")


#### refit rom model with robust variance estimates
secondary.hyp.rom.best.model.prep=rma.mv(yi.rom,vi.rom,data = secondary.data.hyp,
                                       mods = ~ avg.rir+load.set+set.rep.equated+weeks+train.status,
                                       random = list(~1|study,~1|group,~1|obs),
                                       method = "REML", test = "t",dfs = "contain")

secondary.hyp.rom.best.model=robust.rma.mv(secondary.hyp.rom.best.model.prep,
                                         cluster = secondary.data.hyp$study)

#### I2 and R2 for rom model
secondary.hyp.rom.best.model.r2=r2_ml(secondary.hyp.rom.best.model)*100
save(secondary.hyp.rom.best.model.r2,file = "secondary.hyp.rom.best.model.r2.RData")
secondary.hyp.rom.best.model.i2=i2_ml(secondary.hyp.rom.best.model)
save(secondary.hyp.rom.best.model.i2,file = "secondary.hyp.rom.best.model.i2.RData")

#### rom model output
secondary.hyp.rom.best.output<-subset(as.data.frame(coef(summary(secondary.hyp.rom.best.model))),select = c(-pval,-tval))
secondary.hyp.rom.best.output.prep1<-qt(.975,secondary.hyp.rom.best.output$df)
secondary.hyp.rom.best.output.prep2<-sqrt((secondary.hyp.rom.best.output$se^2)+sum(secondary.hyp.rom.best.model$sigma2))
secondary.hyp.rom.best.output$lower.PI<-secondary.hyp.rom.best.output$estimate-(secondary.hyp.rom.best.output.prep1*secondary.hyp.rom.best.output.prep2)
secondary.hyp.rom.best.output$upper.PI<-secondary.hyp.rom.best.output$estimate+(secondary.hyp.rom.best.output.prep1*secondary.hyp.rom.best.output.prep2)
secondary.hyp.rom.best.output$Parameter<-c("Intercept",
                                         "Estimated RIR",
                                         "Load",
                                         "Volume Equated [Rep]",
                                         "Volume Equated [Both]",
                                         "Weeks",
                                         "Training Status [Untrained]")

secondary.hyp.rom.best.output.table<-relocate(secondary.hyp.rom.best.output,Parameter)%>%
  gt()%>%
  cols_label(ci.lb="Lower Bound",
             ci.ub="Upper Bound",
             lower.PI="Lower Bound",
             upper.PI="Upper Bound",
             se="Standard Error",
             estimate="Estimate")%>%
  fmt_markdown(columns = 1)%>%
  tab_header(title = "Coefficients of Linear Multilevel Model for Muscle Hypertrophy with <10 RIR (Response Ratio)")%>%
  tab_spanner(columns = 5:6,label = "95% Confidence Interval")%>%
  tab_spanner(columns = 7:8,label = "95% Predicition Interval")%>%
  opt_stylize(style = 1,color = "gray")%>%
  tab_options(column_labels.background.color ="gray97",column_labels.font.weight = "bold")%>%
  tab_style(style = cell_text(style = "italic"),locations = cells_row_groups())%>%
  tab_style(style = cell_text(align = "center"),locations = cells_column_labels())%>%
  tab_style(style = cell_text(weight  = "bold"),locations = cells_title())%>%
  tab_style(style = cell_text(align = "center"),locations = cells_body(columns=-1))%>%
  cols_width(-df~px(133))%>%
  cols_width(Parameter~px(200))

gtsave(secondary.hyp.rom.best.output.table,filename = "secondary.hyp.rom.best.output.table.png",
       path = "/Users/zacrobinson/Documents/For Publication PTF Meta/Models")

# MAIN EFFECTS ------------------------------------------------------

### strength
#### g marginal means
secondary.str.g.best.emmeans=emmprep(secondary.str.g.best.model,
                                   at=list(avg.rir=0:10))%>%
  emmeans(~avg.rir,
          weights="prop")%>%
  pred_interval_esmeans(model=secondary.str.g.best.model)%>%
  as.data.frame()

#### g marginal slope
secondary.str.g.best.slope.prep<-emmprep(secondary.str.g.best.model,
                                       at=list(avg.rir=secondary.str.mean.rir))

secondary.str.g.best.slope<-confint(emmeans(secondary.str.g.best.slope.prep,consec~avg.rir,
                                          weights = "prop")$contrasts)

secondary.str.g.best.slope.prep2<-qt(.975,secondary.str.g.best.slope$df)
secondary.str.g.best.slope.prep3<-sqrt((secondary.str.g.best.slope$SE^2)+sum(secondary.str.g.best.model$sigma2))

secondary.str.g.best.slope$lower.PI<-secondary.str.g.best.slope$estimate-(secondary.str.g.best.slope.prep2*secondary.str.g.best.slope.prep3)
secondary.str.g.best.slope$upper.PI<-secondary.str.g.best.slope$estimate+(secondary.str.g.best.slope.prep2*secondary.str.g.best.slope.prep3)

save(secondary.str.g.best.slope,file = "secondary.str.g.best.slope.RData")

#### g plot
secondary.str.g.best.model.plot=ggplot(data=secondary.str.g.best.emmeans,
                                     aes(x=avg.rir,
                                         y=emmean))+THEME+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(data = secondary.data.str,aes(x=avg.rir,y=yi,size=weights),alpha=0.75,shape=16,color="#c1272d")+
  geom_ribbon(aes(x=avg.rir,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(aes(x=avg.rir,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(linewidth=1,color="#c1272d")+
  scale_x_continuous(breaks = c(0,2,4,6,8,10))+
  scale_y_continuous(n.breaks = 9)+
  scale_size_continuous(range = c(2,8))+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Marginal Effects from Linear-log Multilevel Model for Maximal Strength with <10 RIR (Standardized Mean Change)")

ggsave(secondary.str.g.best.model.plot,device = "pdf",filename="secondary.str.g.best.model.plot.pdf",width = 12,height = 8,
       path = "/Users/zacrobinson/Documents/For Publication PTF Meta/Plots")

#### g table
secondary.str.g.best.model.table=secondary.str.g.best.emmeans%>%
  gt()%>%
  tab_header(title = "Marginal Means from Linear-log Multilevel Model for Maximal Strength with <10 RIR (Standardized Mean Change)")%>%
  cols_label(avg.rir="Estimated Proximity to Failure (RIR)",
             emmean="Effect Size (Hedge's g)",
             SE="Standard Error",
             lower.CL="Lower Bound",
             upper.CL="Upper Bound",
             lower.PI="Lower Bound",
             upper.PI="Upper Bound")%>%
  tab_spanner(label = "95% Confidence Interval",columns = 5:6)%>%
  tab_spanner(label = "95% Predicition Interval",columns = 7:8)%>%
  opt_stylize(style = 1,color = "gray")%>%
  tab_options(column_labels.background.color ="gray97",column_labels.font.weight = "bold")%>%
  tab_style(style = cell_text(style = "italic"),locations = cells_row_groups())%>%
  tab_style(style = cell_text(align = "center"),locations = cells_column_labels())%>%
  tab_style(style = cell_text(weight  = "bold"),locations = cells_title())%>%
  cols_width(-df~px(150))%>%
  tab_style(style = cell_text(align = "center"),locations = cells_body())

gtsave(secondary.str.g.best.model.table,filename = "secondary.str.g.best.model.table.png",vheight=1000,vwidth=1200,
       path = "/Users/zacrobinson/Documents/For Publication PTF Meta/Tables")

#### rom marginal means
secondary.str.rom.best.emmeans=emmprep(secondary.str.rom.best.model,
                                     at=list(avg.rir=0:10))%>%
  emmeans(~avg.rir,
          weights="prop")%>%
  pred_interval_esmeans(model=secondary.str.rom.best.model)%>%
  as.data.frame()%>%
  mutate_at(c(2,3,5,6,7,8),expo)

#### rom marginal slope
secondary.str.rom.best.slope.prep<-emmprep(secondary.str.rom.best.model,
                                         at=list(avg.rir=secondary.str.mean.rir))

secondary.str.rom.best.slope<-confint(emmeans(secondary.str.rom.best.slope.prep,consec~avg.rir,
                                            weights = "prop")$contrasts)

secondary.str.rom.best.slope.prep2<-qt(.975,secondary.str.rom.best.slope$df)
secondary.str.rom.best.slope.prep3<-sqrt((secondary.str.rom.best.slope$SE^2)+sum(secondary.str.rom.best.model$sigma2))

secondary.str.rom.best.slope$lower.PI<-secondary.str.rom.best.slope$estimate-(secondary.str.rom.best.slope.prep2*secondary.str.rom.best.slope.prep3)
secondary.str.rom.best.slope$upper.PI<-secondary.str.rom.best.slope$estimate+(secondary.str.rom.best.slope.prep2*secondary.str.rom.best.slope.prep3)

secondary.str.rom.best.slope=secondary.str.rom.best.slope%>%
  mutate_at(c(2,3,5,6,7,8),expo)

save(secondary.str.rom.best.slope,file = "secondary.str.rom.best.slope.RData")

#### rom plot
secondary.str.rom.best.model.plot=ggplot(data=secondary.str.rom.best.emmeans,
                                       aes(x=avg.rir,
                                           y=emmean))+THEME+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(data = secondary.data.str,aes(x=avg.rir,y=yi.rom.exp,size=weights.rom),alpha=0.75,shape=16,color="#c1272d")+
  geom_ribbon(aes(x=avg.rir,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(aes(x=avg.rir,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(linewidth=1,color="#c1272d")+
  scale_x_continuous(breaks = c(0,2,4,6,8,10))+
  scale_y_continuous(n.breaks = 9)+
  scale_size_continuous(range = c(2,8))+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio for Muscle Size (% Change)",
       title = "Marginal Effects from Linear Multilevel Model for Maximal Strength with <10 RIR (Response Ratio)")

ggsave(secondary.str.rom.best.model.plot,device = "pdf",filename="secondary.str.rom.best.model.plot.pdf",width = 12,height = 8,
       path = "/Users/zacrobinson/Documents/For Publication PTF Meta/Plots")

#### rom table
secondary.str.rom.best.model.table=secondary.str.rom.best.emmeans%>%
  gt()%>%
  tab_header(title = "Marginal Means from Linear Multilevel Model for Maximal Strength with <10 RIR (Response Ratio)")%>%
  cols_label(avg.rir="Estimated Proximity to Failure (RIR)",
             emmean="Effect Size (% Change)",
             SE="Standard Error",
             lower.CL="Lower Bound",
             upper.CL="Upper Bound",
             lower.PI="Lower Bound",
             upper.PI="Upper Bound")%>%
  tab_spanner(label = "95% Confidence Interval",columns = 5:6)%>%
  tab_spanner(label = "95% Predicition Interval",columns = 7:8)%>%
  opt_stylize(style = 1,color = "gray")%>%
  tab_options(column_labels.background.color ="gray97",column_labels.font.weight = "bold")%>%
  tab_style(style = cell_text(style = "italic"),locations = cells_row_groups())%>%
  tab_style(style = cell_text(align = "center"),locations = cells_column_labels())%>%
  tab_style(style = cell_text(weight  = "bold"),locations = cells_title())%>%
  cols_width(-df~px(150))%>%
  tab_style(style = cell_text(align = "center"),locations = cells_body())

gtsave(secondary.str.rom.best.model.table,filename = "secondary.str.rom.best.model.table.png",vheight=1000,vwidth=1200,
       path = "/Users/zacrobinson/Documents/For Publication PTF Meta/Tables")

### hypertrophy
#### g marginal means
secondary.hyp.g.best.emmeans=emmprep(secondary.hyp.g.best.model,
                                   at=list(avg.rir=0:10))%>%
  emmeans(~avg.rir,
          weights="prop")%>%
  pred_interval_esmeans(model=secondary.hyp.g.best.model)%>%
  as.data.frame()

#### g marginal slope
secondary.hyp.g.best.slope.prep<-emmprep(secondary.hyp.g.best.model,
                                       at=list(avg.rir=secondary.hyp.mean.rir))

secondary.hyp.g.best.slope<-confint(emmeans(secondary.hyp.g.best.slope.prep,consec~avg.rir,
                                          weights = "prop")$contrasts)

secondary.hyp.g.best.slope.prep2<-qt(.975,secondary.hyp.g.best.slope$df)
secondary.hyp.g.best.slope.prep3<-sqrt((secondary.hyp.g.best.slope$SE^2)+sum(secondary.hyp.g.best.model$sigma2))

secondary.hyp.g.best.slope$lower.PI<-secondary.hyp.g.best.slope$estimate-(secondary.hyp.g.best.slope.prep2*secondary.hyp.g.best.slope.prep3)
secondary.hyp.g.best.slope$upper.PI<-secondary.hyp.g.best.slope$estimate+(secondary.hyp.g.best.slope.prep2*secondary.hyp.g.best.slope.prep3)

save(secondary.hyp.g.best.slope,file = "secondary.hyp.g.best.slope.RData")

#### g plot
secondary.hyp.g.best.model.plot=ggplot(data=secondary.hyp.g.best.emmeans,
                                     aes(x=avg.rir,
                                         y=emmean))+THEME+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(data = secondary.data.hyp,aes(x=avg.rir,y=yi,size=weights),alpha=0.75,shape=16,color="#0000a7")+
  geom_ribbon(aes(x=avg.rir,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(aes(x=avg.rir,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(linewidth=1,color="#0000a7")+
  scale_x_continuous(breaks = c(0,2,4,6,8,10))+
  scale_y_continuous(n.breaks = 9)+
  scale_size_continuous(range = c(2,8))+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Muscle Size (Hedge's g)",
       title = "Marginal Effects from Linear Multilevel Model for Muscle Hypertrophy with <10 RIR (Standardized Mean Change)")

ggsave(secondary.hyp.g.best.model.plot,device = "pdf",filename="secondary.hyp.g.best.model.plot.pdf",width = 12,height = 8,
       path = "/Users/zacrobinson/Documents/For Publication PTF Meta/Plots")

#### g table
secondary.hyp.g.best.model.table=secondary.hyp.g.best.emmeans%>%
  gt()%>%
  tab_header(title = "Marginal Means from Linear Multilevel Model for Muscle Hypertrophy with <10 RIR (Standardized Mean Change)")%>%
  cols_label(avg.rir="Estimated Proximity to Failure (RIR)",
             emmean="Effect Size (Hedge's g)",
             SE="Standard Error",
             lower.CL="Lower Bound",
             upper.CL="Upper Bound",
             lower.PI="Lower Bound",
             upper.PI="Upper Bound")%>%
  tab_spanner(label = "95% Confidence Interval",columns = 5:6)%>%
  tab_spanner(label = "95% Predicition Interval",columns = 7:8)%>%
  opt_stylize(style = 1,color = "gray")%>%
  tab_options(column_labels.background.color ="gray97",column_labels.font.weight = "bold")%>%
  tab_style(style = cell_text(style = "italic"),locations = cells_row_groups())%>%
  tab_style(style = cell_text(align = "center"),locations = cells_column_labels())%>%
  tab_style(style = cell_text(weight  = "bold"),locations = cells_title())%>%
  cols_width(-df~px(150))%>%
  tab_style(style = cell_text(align = "center"),locations = cells_body())

gtsave(secondary.hyp.g.best.model.table,filename = "secondary.hyp.g.best.model.table.png",vheight=1000,vwidth=1200,
       path = "/Users/zacrobinson/Documents/For Publication PTF Meta/Tables")



#### rom marginal means
secondary.hyp.rom.best.emmeans=emmprep(secondary.hyp.rom.best.model,
                                     at=list(avg.rir=0:10))%>%
  emmeans(~avg.rir,
          weights="prop")%>%
  pred_interval_esmeans(model=secondary.hyp.rom.best.model)%>%
  as.data.frame()%>%
  mutate_at(c(2,3,5,6,7,8),expo)

#### rom marginal slope
secondary.hyp.rom.best.slope.prep<-emmprep(secondary.hyp.rom.best.model,
                                         at=list(avg.rir=secondary.hyp.mean.rir))

secondary.hyp.rom.best.slope<-confint(emmeans(secondary.hyp.rom.best.slope.prep,consec~avg.rir,
                                            weights = "prop")$contrasts)

secondary.hyp.rom.best.slope.prep2<-qt(.975,secondary.hyp.rom.best.slope$df)
secondary.hyp.rom.best.slope.prep3<-sqrt((secondary.hyp.rom.best.slope$SE^2)+sum(secondary.hyp.rom.best.model$sigma2))

secondary.hyp.rom.best.slope$lower.PI<-secondary.hyp.rom.best.slope$estimate-(secondary.hyp.rom.best.slope.prep2*secondary.hyp.rom.best.slope.prep3)
secondary.hyp.rom.best.slope$upper.PI<-secondary.hyp.rom.best.slope$estimate+(secondary.hyp.rom.best.slope.prep2*secondary.hyp.rom.best.slope.prep3)

secondary.hyp.rom.best.slope=secondary.hyp.rom.best.slope%>%
  mutate_at(c(2,3,5,6,7,8),expo)

save(secondary.hyp.rom.best.slope,file = "secondary.hyp.rom.best.slope.RData")

#### rom plot
secondary.hyp.rom.best.model.plot=ggplot(data=secondary.hyp.rom.best.emmeans,
                                       aes(x=avg.rir,
                                           y=emmean))+THEME+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(data = secondary.data.hyp,aes(x=avg.rir,y=yi.rom.exp,size=weights.rom),alpha=0.75,shape=16,color="#0000a7")+
  geom_ribbon(aes(x=avg.rir,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(aes(x=avg.rir,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(linewidth=1,color="#0000a7")+
  scale_x_continuous(breaks = c(0,2,4,6,8,10))+
  scale_y_continuous(n.breaks = 9)+
  scale_size_continuous(range = c(2,8))+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio for Muscle Size (% Change)",
       title = "Marginal Effects from Linear Multilevel Model for Muscle Hypertrophy with <10 RIR (Response Ratio)")

ggsave(secondary.hyp.rom.best.model.plot,device = "pdf",filename="secondary.hyp.rom.best.model.plot.pdf",width = 12,height = 8,
       path = "/Users/zacrobinson/Documents/For Publication PTF Meta/Plots")

#### rom table
secondary.hyp.rom.best.model.table=secondary.hyp.rom.best.emmeans%>%
  gt()%>%
  tab_header(title = "Marginal Means from Linear Multilevel Model for Muscle Hypertrophy with <10 RIR (Response Ratio)")%>%
  cols_label(avg.rir="Estimated Proximity to Failure (RIR)",
             emmean="Effect Size (% Change)",
             SE="Standard Error",
             lower.CL="Lower Bound",
             upper.CL="Upper Bound",
             lower.PI="Lower Bound",
             upper.PI="Upper Bound")%>%
  tab_spanner(label = "95% Confidence Interval",columns = 5:6)%>%
  tab_spanner(label = "95% Predicition Interval",columns = 7:8)%>%
  opt_stylize(style = 1,color = "gray")%>%
  tab_options(column_labels.background.color ="gray97",column_labels.font.weight = "bold")%>%
  tab_style(style = cell_text(style = "italic"),locations = cells_row_groups())%>%
  tab_style(style = cell_text(align = "center"),locations = cells_column_labels())%>%
  tab_style(style = cell_text(weight  = "bold"),locations = cells_title())%>%
  cols_width(-df~px(150))%>%
  tab_style(style = cell_text(align = "center"),locations = cells_body())

gtsave(secondary.hyp.rom.best.model.table,filename = "secondary.hyp.rom.best.model.table.png",vheight=1000,vwidth=1200,
       path = "/Users/zacrobinson/Documents/For Publication PTF Meta/Tables")

# MODERATORS --------------------------------------------------------------

### strength
#### g load
secondary.str.g.moderator.load=secondary.str.g.best.model.prep%>%
  update(~. + log1p(avg.rir):load.set)%>%
  robust.rma.mv(cluster = secondary.data.str$study)
  
secondary.str.g.moderator.load.emms=secondary.str.g.moderator.load%>%
  emmprep(at=list(avg.rir=0:10,
                  load.set=c(30,60,90)))%>%
  emmeans(~avg.rir+load.set,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.str.g.moderator.load)

secondary.str.g.moderator.load.contrast=emmprep(secondary.str.g.moderator.load,
                                              at=list(avg.rir=secondary.str.mean.rir,
                                                      load.set=secondary.str.mean.load.set))%>%
  emmeans(~avg.rir+load.set,
          weights = "prop")%>%
  contrast("consec",interaction=TRUE)%>%
  confint()

secondary.str.g.moderator.load.contrast.label=paste0('Interaction Contrast (Load) = ',
                                                   bquote(.(round(secondary.str.g.moderator.load.contrast[[3]],5))),
                                                   " [95% CI: ",bquote(.(round(secondary.str.g.moderator.load.contrast[[6]],5))),
                                                   ", ",bquote(.(round(secondary.str.g.moderator.load.contrast[[7]],5))),"]")
  
secondary.str.g.moderator.load.plot=ggplot(data = secondary.data.str%>%
         mutate(load.set=ifelse(load.set<45,"30% of 1RM",
                                ifelse(load.set>45&load.set<75,"60% of 1RM","90% of 1RM"))),
       aes(x=avg.rir,
           y=yi))+THEME+theme(plot.caption = element_text(face = "italic"))+
  facet_wrap(~load.set)+
  labs(x="Estimated Proximity to Failure (RIR)",
                            y="Standardized Mean Change in Strength (Hedge's g)",
                            title = "Load (% of 1RM) x <10 Estimated RIR",
       caption = secondary.str.g.moderator.load.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = secondary.str.g.moderator.load.emms%>%
                mutate(load.set=ifelse(load.set==30,"30% of 1RM",
                                       ifelse(load.set==60,"60% of 1RM","90% of 1RM"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = secondary.str.g.moderator.load.emms%>%
                mutate(load.set=ifelse(load.set==30,"30% of 1RM",
                                       ifelse(load.set==60,"60% of 1RM","90% of 1RM"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = secondary.str.g.moderator.load.emms%>%
              mutate(load.set=ifelse(load.set==30,"30% of 1RM",
                                     ifelse(load.set==60,"60% of 1RM","90% of 1RM"))),
    aes(x=avg.rir,
        y=emmean),
    linewidth=1,
    color="#c1272d")+
  scale_size_continuous(guide="none")

#### g volume equating method
secondary.str.g.moderator.volume.method=secondary.str.g.best.model.prep%>%
  update(~. + log1p(avg.rir):set.rep.equated)%>%
  robust.rma.mv(cluster = secondary.data.str$study)

secondary.str.g.moderator.volume.method.emms=secondary.str.g.moderator.volume.method%>%emmprep(at=list(avg.rir=0:10))%>%
  emmeans(~avg.rir+set.rep.equated,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.str.g.moderator.volume.method)

secondary.str.g.moderator.volume.method.contrast=emmprep(secondary.str.g.moderator.volume.method,
                                              at=list(avg.rir=secondary.str.mean.rir))%>%
  emmeans(~avg.rir+set.rep.equated,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()

secondary.str.g.moderator.volume.method.contrast.label=paste0('Interaction Contrast (Set - Rep) = ',
                                                   bquote(.(round(secondary.str.g.moderator.volume.method.contrast[[1,3]],5))),
                                                   " [95% CI: ",bquote(.(round(secondary.str.g.moderator.volume.method.contrast[[1,6]],5))),
                                                   ", ",bquote(.(round(secondary.str.g.moderator.volume.method.contrast[[1,7]],5))),"]\n",
                                                   'Interaction Contrast (Set - Both) = ',
                                                   bquote(.(round(secondary.str.g.moderator.volume.method.contrast[[2,3]],5))),
                                                   " [95% CI: ",bquote(.(round(secondary.str.g.moderator.volume.method.contrast[[2,6]],5))),
                                                   ", ",bquote(.(round(secondary.str.g.moderator.volume.method.contrast[[2,7]],5))),"]\n",
                                                   'Interaction Contrast (Rep - Both) = ',
                                                   bquote(.(round(secondary.str.g.moderator.volume.method.contrast[[3,3]],5))),
                                                   " [95% CI: ",bquote(.(round(secondary.str.g.moderator.volume.method.contrast[[3,6]],5))),
                                                   ", ",bquote(.(round(secondary.str.g.moderator.volume.method.contrast[[3,7]],5))),"]")


secondary.str.g.moderator.volume.method.plot=ggplot(data = secondary.data.str%>%
         mutate(set.rep.equated=ifelse(set.rep.equated=="set","Set Equated",
                                       ifelse(set.rep.equated=="rep","Repetition Equated","Both Equated"))),
       aes(x=avg.rir,
           y=yi))+THEME+theme(legend.position = "right",
                              plot.caption = element_text(face = "italic"))+
  facet_wrap(~set.rep.equated)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Method of Volume Equating x <10 Estimated RIR",
       caption = secondary.str.g.moderator.volume.method.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = secondary.str.g.moderator.volume.method.emms%>%
                mutate(set.rep.equated=ifelse(set.rep.equated=="set","Set Equated",
                                              ifelse(set.rep.equated=="rep","Repetition Equated","Both Equated"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = secondary.str.g.moderator.volume.method.emms%>%
                mutate(set.rep.equated=ifelse(set.rep.equated=="set","Set Equated",
                                              ifelse(set.rep.equated=="rep","Repetition Equated","Both Equated"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = secondary.str.g.moderator.volume.method.emms%>%
              mutate(set.rep.equated=ifelse(set.rep.equated=="set","Set Equated",
                                            ifelse(set.rep.equated=="rep","Repetition Equated","Both Equated"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### g training status
secondary.str.g.moderator.train.status=secondary.str.g.best.model.prep%>%
  update(~. + log1p(avg.rir):train.status)%>%
  robust.rma.mv(cluster = secondary.data.str$study)

secondary.str.g.moderator.train.status.emms=secondary.str.g.moderator.train.status%>%emmprep(at=list(avg.rir=0:10))%>%
  emmeans(~avg.rir+train.status,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.str.g.moderator.train.status)

secondary.str.g.moderator.train.status.contrast=emmprep(secondary.str.g.moderator.train.status,
                                                       at=list(avg.rir=secondary.str.mean.rir))%>%
  emmeans(~avg.rir+train.status,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()

secondary.str.g.moderator.train.status.contrast.label=paste0('Interaction Contrast (Trained - Untrained) = ',
                                                            bquote(.(round(secondary.str.g.moderator.train.status.contrast[[1,3]],5))),
                                                            " [95% CI: ",bquote(.(round(secondary.str.g.moderator.train.status.contrast[[1,6]],5))),
                                                            ", ",bquote(.(round(secondary.str.g.moderator.train.status.contrast[[1,7]],5))),"]")


secondary.str.g.moderator.train.status.plot=ggplot(data = secondary.data.str%>%
                                                    mutate(train.status=ifelse(train.status=="trained","Trained","Untrained")),
                                                  aes(x=avg.rir,
                                                      y=yi))+THEME+theme(legend.position = "right",
                                                                         plot.caption = element_text(face = "italic"))+
  facet_wrap(~train.status)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Training Status x <10 Estimated RIR",
       caption = secondary.str.g.moderator.train.status.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = secondary.str.g.moderator.train.status.emms%>%
                mutate(train.status=ifelse(train.status=="trained","Trained","Untrained")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = secondary.str.g.moderator.train.status.emms%>%
                mutate(train.status=ifelse(train.status=="trained","Trained","Untrained")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = secondary.str.g.moderator.train.status.emms%>%
              mutate(train.status=ifelse(train.status=="trained","Trained","Untrained")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### g weeks
secondary.str.g.moderator.weeks=secondary.str.g.best.model.prep%>%
  update(~. + log1p(avg.rir):weeks)%>%
  robust.rma.mv(cluster = secondary.data.str$study)

secondary.str.g.moderator.weeks.emms=secondary.str.g.moderator.weeks%>%
  emmprep(at=list(avg.rir=0:10,
                  weeks=c(4,8,12)))%>%
  emmeans(~avg.rir+weeks,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.str.g.moderator.weeks)

secondary.str.g.moderator.weeks.contrast=emmprep(secondary.str.g.moderator.weeks,
                                              at=list(avg.rir=secondary.str.mean.rir,
                                                      weeks=secondary.str.mean.weeks))%>%
  emmeans(~avg.rir+weeks,
          weights = "prop")%>%
  contrast("consec",interaction=TRUE)%>%
  confint()

secondary.str.g.moderator.weeks.contrast.label=paste0('Interaction Contrast (Weeks) = ',
                                                   bquote(.(round(secondary.str.g.moderator.weeks.contrast[[3]],5))),
                                                   " [95% CI: ",bquote(.(round(secondary.str.g.moderator.weeks.contrast[[6]],5))),
                                                   ", ",bquote(.(round(secondary.str.g.moderator.weeks.contrast[[7]],5))),"]")

secondary.str.g.moderator.weeks.plot=ggplot(data = secondary.data.str%>%
                                           mutate(weeks=ifelse(weeks<6,"4 Weeks",
                                                                  ifelse(weeks>6&weeks<10,"8 Weeks","12 Weeks"))),
                                         aes(x=avg.rir,
                                             y=yi))+THEME+theme(plot.caption = element_text(face = "italic"))+
  facet_wrap(~fct_relevel(weeks,c("4 Weeks","8 Weeks","12 Weeks")))+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Intervention Duration (Weeks) x <10 Estimated RIR",
       caption = secondary.str.g.moderator.weeks.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = secondary.str.g.moderator.weeks.emms%>%
                mutate(weeks=ifelse(weeks==4,"4 Weeks",
                                       ifelse(weeks==8,"8 Weeks","12 Weeks"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = secondary.str.g.moderator.weeks.emms%>%
                mutate(weeks=ifelse(weeks==4,"4 Weeks",
                                    ifelse(weeks==8,"8 Weeks","12 Weeks"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = secondary.str.g.moderator.weeks.emms%>%
              mutate(weeks=ifelse(weeks==4,"4 Weeks",
                                  ifelse(weeks==8,"8 Weeks","12 Weeks"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### g adjusted sets 
secondary.str.g.moderator.adj.sets=secondary.str.g.best.model.prep%>%
  update(~. + log1p(avg.rir)*adj.sets.week)%>%
  robust.rma.mv(cluster = secondary.data.str$study)

secondary.str.g.moderator.adj.sets.emms=secondary.str.g.moderator.adj.sets%>%
  emmprep(at=list(avg.rir=0:10,
                  adj.sets.week=c(8,16,24)))%>%
  emmeans(~avg.rir+adj.sets.week,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.str.g.moderator.adj.sets)

secondary.str.g.moderator.adj.sets.contrast=emmprep(secondary.str.g.moderator.adj.sets,
                                               at=list(avg.rir=secondary.str.mean.rir,
                                                       adj.sets.week=secondary.str.mean.adj.sets.week))%>%
  emmeans(~avg.rir+adj.sets.week,
          weights = "prop")%>%
  contrast("consec",interaction=TRUE)%>%
  confint()

secondary.str.g.moderator.adj.sets.contrast.label=paste0('Interaction Contrast (Adjusted Direct Sets) = ',
                                                    bquote(.(round(secondary.str.g.moderator.adj.sets.contrast[[3]],5))),
                                                    " [95% CI: ",bquote(.(round(secondary.str.g.moderator.adj.sets.contrast[[6]],5))),
                                                    ", ",bquote(.(round(secondary.str.g.moderator.adj.sets.contrast[[7]],5))),"]")

secondary.str.g.moderator.adj.sets.plot=ggplot(data = secondary.data.str%>%
                                            mutate(adj.sets.week=ifelse(adj.sets.week<12,"8 Sets",
                                                                ifelse(adj.sets.week>12&adj.sets.week<20,"16 Sets","24 Sets"))),
                                          aes(x=avg.rir,
                                              y=yi))+THEME+theme(plot.caption = element_text(face = "italic"))+
  facet_wrap(~fct_relevel(adj.sets.week,c("8 Sets","16 Sets","24 Sets")))+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Adjusted Direct Sets Per Week Per Muscle Group x <10 Estimated RIR",
       caption = secondary.str.g.moderator.adj.sets.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = secondary.str.g.moderator.adj.sets.emms%>%
                mutate(adj.sets.week=ifelse(adj.sets.week==8,"8 Sets",
                                       ifelse(adj.sets.week==16,"16 Sets","24 Sets"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = secondary.str.g.moderator.adj.sets.emms%>%
                mutate(adj.sets.week=ifelse(adj.sets.week==8,"8 Sets",
                                       ifelse(adj.sets.week==16,"16 Sets","24 Sets"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = secondary.str.g.moderator.adj.sets.emms%>%
              mutate(adj.sets.week=ifelse(adj.sets.week==8,"8 Sets",
                                     ifelse(adj.sets.week==16,"16 Sets","24 Sets"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### g frequency
secondary.str.g.moderator.frequency.per.muscle=secondary.str.g.best.model.prep%>%
  update(~. + log1p(avg.rir)*frequency.per.muscle)%>%
  robust.rma.mv(cluster = secondary.data.str$study)

secondary.str.g.moderator.frequency.per.muscle.emms=secondary.str.g.moderator.frequency.per.muscle%>%
  emmprep(at=list(avg.rir=0:10,
                  frequency.per.muscle=c(1,2,3)))%>%
  emmeans(~avg.rir+frequency.per.muscle,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.str.g.moderator.frequency.per.muscle)

secondary.str.g.moderator.frequency.per.muscle.contrast=emmprep(secondary.str.g.moderator.frequency.per.muscle,
                                                  at=list(avg.rir=secondary.str.mean.rir,
                                                          frequency.per.muscle=secondary.str.mean.frequency.per.muscle))%>%
  emmeans(~avg.rir+frequency.per.muscle,
          weights = "prop")%>%
  contrast("consec",interaction=TRUE)%>%
  confint()

secondary.str.g.moderator.frequency.per.muscle.contrast.label=paste0('Interaction Contrast (Direct Frequency) = ',
                                                       bquote(.(round(secondary.str.g.moderator.frequency.per.muscle.contrast[[3]],5))),
                                                       " [95% CI: ",bquote(.(round(secondary.str.g.moderator.frequency.per.muscle.contrast[[6]],5))),
                                                       ", ",bquote(.(round(secondary.str.g.moderator.frequency.per.muscle.contrast[[7]],5))),"]")

secondary.str.g.moderator.frequency.per.muscle.plot=ggplot(data = secondary.data.str%>%
                                               mutate(frequency.per.muscle=ifelse(frequency.per.muscle<1.5,"1x Per Week",
                                                                           ifelse(frequency.per.muscle>1.5&frequency.per.muscle<2.5,"2x Per Week","3x Per Week"))),
                                             aes(x=avg.rir,
                                                 y=yi))+THEME+theme(plot.caption = element_text(face = "italic"))+
  facet_wrap(~frequency.per.muscle)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Direct Frequency Per Week Per Muscle Group x <10 Estimated RIR",
       caption = secondary.str.g.moderator.frequency.per.muscle.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = secondary.str.g.moderator.frequency.per.muscle.emms%>%
                mutate(frequency.per.muscle=ifelse(frequency.per.muscle==1,"1x Per Week",
                                       ifelse(frequency.per.muscle==2,"2x Per Week","3x Per Week"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = secondary.str.g.moderator.frequency.per.muscle.emms%>%
                mutate(frequency.per.muscle=ifelse(frequency.per.muscle==1,"1x Per Week",
                                       ifelse(frequency.per.muscle==2,"2x Per Week","3x Per Week"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = secondary.str.g.moderator.frequency.per.muscle.emms%>%
              mutate(frequency.per.muscle=ifelse(frequency.per.muscle==1,"1x Per Week",
                                     ifelse(frequency.per.muscle==2,"2x Per Week","3x Per Week"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### g age
secondary.str.g.moderator.age=secondary.str.g.best.model.prep%>%
  update(~. + log1p(avg.rir)*age)%>%
  robust.rma.mv(cluster = secondary.data.str$study)

secondary.str.g.moderator.age.emms=secondary.str.g.moderator.age%>%
  emmprep(at=list(avg.rir=0:10,
                  age=c(20,40,60)))%>%
  emmeans(~avg.rir+age,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.str.g.moderator.age)

secondary.str.g.moderator.age.contrast=emmprep(secondary.str.g.moderator.age,
                                                              at=list(avg.rir=secondary.str.mean.rir,
                                                                      age=secondary.str.mean.age))%>%
  emmeans(~avg.rir+age,
          weights = "prop")%>%
  contrast("consec",interaction=TRUE)%>%
  confint()

secondary.str.g.moderator.age.contrast.label=paste0('Interaction Contrast (Age) = ',
                                                                   bquote(.(round(secondary.str.g.moderator.age.contrast[[3]],5))),
                                                                   " [95% CI: ",bquote(.(round(secondary.str.g.moderator.age.contrast[[6]],5))),
                                                                   ", ",bquote(.(round(secondary.str.g.moderator.age.contrast[[7]],5))),"]")

secondary.str.g.moderator.age.plot=ggplot(data = secondary.data.str%>%
                                                           mutate(age=ifelse(age<30,"20 Years Old",
                                                                                              ifelse(age>30&age<50,"40 Years Old","60 Years Old"))),
                                                         aes(x=avg.rir,
                                                             y=yi))+THEME+theme(plot.caption = element_text(face = "italic"))+
  facet_wrap(~age)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Age x <10 Estimated RIR",
       caption = secondary.str.g.moderator.age.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = secondary.str.g.moderator.age.emms%>%
                mutate(age=ifelse(age==20,"20 Years Old",
                                       ifelse(age==40,"40 Years Old","60 Years Old"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = secondary.str.g.moderator.age.emms%>%
                mutate(age=ifelse(age==20,"20 Years Old",
                                       ifelse(age==40,"40 Years Old","60 Years Old"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = secondary.str.g.moderator.age.emms%>%
              mutate(age=ifelse(age==20,"20 Years Old",
                                     ifelse(age==40,"40 Years Old","60 Years Old"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### g sex
secondary.str.g.moderator.sex.percent.male=secondary.str.g.best.model.prep%>%
  update(~. + log1p(avg.rir)*sex.percent.male)%>%
  robust.rma.mv(cluster = secondary.data.str$study)

secondary.str.g.moderator.sex.percent.male.emms=secondary.str.g.moderator.sex.percent.male%>%
  emmprep(at=list(avg.rir=0:10,
                  sex.percent.male=c(0,100)))%>%
  emmeans(~avg.rir+sex.percent.male,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.str.g.moderator.sex.percent.male)

secondary.str.g.moderator.sex.percent.male.contrast=emmprep(secondary.str.g.moderator.sex.percent.male,
                                             at=list(avg.rir=secondary.str.mean.rir,
                                                     sex.percent.male=secondary.str.mean.sex))%>%
  emmeans(~avg.rir+sex.percent.male,
          weights = "prop")%>%
  contrast("consec",interaction=TRUE)%>%
  confint()

secondary.str.g.moderator.sex.percent.male.contrast.label=paste0('Interaction Contrast (Sex) = ',
                                                  bquote(.(round(secondary.str.g.moderator.sex.percent.male.contrast[[3]],5))),
                                                  " [95% CI: ",bquote(.(round(secondary.str.g.moderator.sex.percent.male.contrast[[6]],5))),
                                                  ", ",bquote(.(round(secondary.str.g.moderator.sex.percent.male.contrast[[7]],5))),"]")

secondary.str.g.moderator.sex.percent.male.plot=ggplot(data = secondary.data.str%>%
                                          mutate(sex.percent.male=ifelse(sex.percent.male<50,"Female","Male")),
                                        aes(x=avg.rir,
                                            y=yi))+THEME+theme(plot.caption = element_text(face = "italic"))+
  facet_wrap(~sex.percent.male)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Sex (% Male) x <10 Estimated RIR",
       caption = secondary.str.g.moderator.sex.percent.male.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = secondary.str.g.moderator.sex.percent.male.emms%>%
                mutate(sex.percent.male=ifelse(sex.percent.male==0,"Female","Male")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = secondary.str.g.moderator.sex.percent.male.emms%>%
                mutate(sex.percent.male=ifelse(sex.percent.male==0,"Female","Male")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = secondary.str.g.moderator.sex.percent.male.emms%>%
              mutate(sex.percent.male=ifelse(sex.percent.male==0,"Female","Male")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### g study design
secondary.str.g.moderator.within.between.design=secondary.str.g.best.model.prep%>%
  update(~. + log1p(avg.rir)*within.between.design)%>%
  robust.rma.mv(cluster = secondary.data.str$study)

secondary.str.g.moderator.within.between.design.emms=secondary.str.g.moderator.within.between.design%>%emmprep(at=list(avg.rir=0:10))%>%
  emmeans(~avg.rir+within.between.design,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.str.g.moderator.within.between.design)

secondary.str.g.moderator.within.between.design.contrast=emmprep(secondary.str.g.moderator.within.between.design,
                                                      at=list(avg.rir=secondary.str.mean.rir))%>%
  emmeans(~avg.rir+within.between.design,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()

secondary.str.g.moderator.within.between.design.contrast.label=paste0('Interaction Contrast (Between - Within) = ',
                                                           bquote(.(round(secondary.str.g.moderator.within.between.design.contrast[[1,3]],5))),
                                                           " [95% CI: ",bquote(.(round(secondary.str.g.moderator.within.between.design.contrast[[1,6]],5))),
                                                           ", ",bquote(.(round(secondary.str.g.moderator.within.between.design.contrast[[1,7]],5))),"]")


secondary.str.g.moderator.within.between.design.plot=ggplot(data = secondary.data.str%>%
                                                   mutate(within.between.design=ifelse(within.between.design=="Between","Between Subjects","Within Subjects")),
                                                 aes(x=avg.rir,
                                                     y=yi))+THEME+theme(legend.position = "right",
                                                                        plot.caption = element_text(face = "italic"))+
  facet_wrap(~within.between.design)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Study Design x <10 Estimated RIR",
       caption = secondary.str.g.moderator.within.between.design.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = secondary.str.g.moderator.within.between.design.emms%>%
                mutate(within.between.design=ifelse(within.between.design=="Between","Between Subjects","Within Subjects")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = secondary.str.g.moderator.within.between.design.emms%>%
                mutate(within.between.design=ifelse(within.between.design=="Between","Between Subjects","Within Subjects")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = secondary.str.g.moderator.within.between.design.emms%>%
              mutate(within.between.design=ifelse(within.between.design=="Between","Between Subjects","Within Subjects")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### g upper / lower outcome
secondary.str.g.moderator.upper.lower.other=secondary.str.g.best.model.prep%>%
  update(~. + log1p(avg.rir)*upper.lower.other)%>%
  robust.rma.mv(cluster = secondary.data.str$study)

secondary.str.g.moderator.upper.lower.other.emms=secondary.str.g.moderator.upper.lower.other%>%emmprep(at=list(avg.rir=0:10))%>%
  emmeans(~avg.rir+upper.lower.other,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.str.g.moderator.upper.lower.other)

secondary.str.g.moderator.upper.lower.other.contrast=emmprep(secondary.str.g.moderator.upper.lower.other,
                                                               at=list(avg.rir=secondary.str.mean.rir))%>%
  emmeans(~avg.rir+upper.lower.other,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()

secondary.str.g.moderator.upper.lower.other.contrast.label=paste0('Interaction Contrast (Lower - Upper) = ',
                                                                    bquote(.(round(secondary.str.g.moderator.upper.lower.other.contrast[[1,3]],5))),
                                                                    " [95% CI: ",bquote(.(round(secondary.str.g.moderator.upper.lower.other.contrast[[1,6]],5))),
                                                                    ", ",bquote(.(round(secondary.str.g.moderator.upper.lower.other.contrast[[1,7]],5))),"]")


secondary.str.g.moderator.upper.lower.other.plot=ggplot(data = secondary.data.str%>%
                                                            mutate(upper.lower.other=ifelse(upper.lower.other=="Lower","Lower Body Outcomes","Upper Body Outcomes")),
                                                          aes(x=avg.rir,
                                                              y=yi))+THEME+theme(legend.position = "right",
                                                                                 plot.caption = element_text(face = "italic"))+
  facet_wrap(~upper.lower.other)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Upper/Lower Body Outcome x <10 Estimated RIR",
       caption = secondary.str.g.moderator.upper.lower.other.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = secondary.str.g.moderator.upper.lower.other.emms%>%
                mutate(upper.lower.other=ifelse(upper.lower.other=="Lower","Lower Body Outcomes","Upper Body Outcomes")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = secondary.str.g.moderator.upper.lower.other.emms%>%
                mutate(upper.lower.other=ifelse(upper.lower.other=="Lower","Lower Body Outcomes","Upper Body Outcomes")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = secondary.str.g.moderator.upper.lower.other.emms%>%
              mutate(upper.lower.other=ifelse(upper.lower.other=="Lower","Lower Body Outcomes","Upper Body Outcomes")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### g exercise selection
secondary.str.g.moderator.train.exercise=secondary.str.g.best.model.prep%>%
  update(~. + log1p(avg.rir)*train.exercise)%>%
  robust.rma.mv(cluster = secondary.data.str$study)

secondary.str.g.moderator.train.exercise.emms=secondary.str.g.moderator.train.exercise%>%emmprep(at=list(avg.rir=0:10))%>%
  emmeans(~avg.rir+train.exercise,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.str.g.moderator.train.exercise)

secondary.str.g.moderator.train.exercise.contrast=emmprep(secondary.str.g.moderator.train.exercise,
                                                           at=list(avg.rir=secondary.str.mean.rir))%>%
  emmeans(~avg.rir+train.exercise,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()

secondary.str.g.moderator.train.exercise.contrast.label=paste0('Interaction Contrast (Both - Multi-joint) = ',
                                                                bquote(.(round(secondary.str.g.moderator.train.exercise.contrast[[1,3]],5))),
                                                                " [95% CI: ",bquote(.(round(secondary.str.g.moderator.train.exercise.contrast[[1,6]],5))),
                                                                ", ",bquote(.(round(secondary.str.g.moderator.train.exercise.contrast[[1,7]],5))),"]\n",
                                                                'Interaction Contrast (Both - Single-joint) = ',
                                                                bquote(.(round(secondary.str.g.moderator.train.exercise.contrast[[2,3]],5))),
                                                                " [95% CI: ",bquote(.(round(secondary.str.g.moderator.train.exercise.contrast[[2,6]],5))),
                                                                ", ",bquote(.(round(secondary.str.g.moderator.train.exercise.contrast[[2,7]],5))),"]\n",
                                                             'Interaction Contrast (Multi-joint - Single-joint) = ',
                                                             bquote(.(round(secondary.str.g.moderator.train.exercise.contrast[[3,3]],5))),
                                                             " [95% CI: ",bquote(.(round(secondary.str.g.moderator.train.exercise.contrast[[3,6]],5))),
                                                             ", ",bquote(.(round(secondary.str.g.moderator.train.exercise.contrast[[3,7]],5))),"]")


secondary.str.g.moderator.train.exercise.plot=ggplot(data = secondary.data.str%>%
                                                        mutate(train.exercise=ifelse(train.exercise=="Both","Both",
                                                                                     ifelse(train.exercise=="Multi-joint","Multi-joint Exercises","Single-joint Exercises"))),
                                                      aes(x=avg.rir,
                                                          y=yi))+THEME+theme(legend.position = "right",
                                                                             plot.caption = element_text(face = "italic"))+
  facet_wrap(~train.exercise)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Exercise Selection x <10 Estimated RIR",
       caption = secondary.str.g.moderator.train.exercise.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = secondary.str.g.moderator.train.exercise.emms%>%
                mutate(train.exercise=ifelse(train.exercise=="Both","Both",
                                             ifelse(train.exercise=="Multi-joint","Multi-joint Exercises","Single-joint Exercises"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = secondary.str.g.moderator.train.exercise.emms%>%
                mutate(train.exercise=ifelse(train.exercise=="Both","Both",
                                             ifelse(train.exercise=="Multi-joint","Multi-joint Exercises","Single-joint Exercises"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = secondary.str.g.moderator.train.exercise.emms%>%
              mutate(train.exercise=ifelse(train.exercise=="Both","Both",
                                           ifelse(train.exercise=="Multi-joint","Multi-joint Exercises","Single-joint Exercises"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### g concurrent training
secondary.str.g.moderator.formal.cardio.intervention=secondary.str.g.best.model.prep%>%
  update(~. + log1p(avg.rir)*formal.cardio.intervention)%>%
  robust.rma.mv(cluster = secondary.data.str$study)

secondary.str.g.moderator.formal.cardio.intervention.emms=secondary.str.g.moderator.formal.cardio.intervention%>%emmprep(at=list(avg.rir=0:10))%>%
  emmeans(~avg.rir+formal.cardio.intervention,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.str.g.moderator.formal.cardio.intervention)

secondary.str.g.moderator.formal.cardio.intervention.contrast=emmprep(secondary.str.g.moderator.formal.cardio.intervention,
                                                           at=list(avg.rir=secondary.str.mean.rir))%>%
  emmeans(~avg.rir+formal.cardio.intervention,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()

secondary.str.g.moderator.formal.cardio.intervention.contrast.label=paste0('Interaction Contrast (RT Only - Concurrent) = ',
                                                                bquote(.(round(secondary.str.g.moderator.formal.cardio.intervention.contrast[[1,3]],5))),
                                                                " [95% CI: ",bquote(.(round(secondary.str.g.moderator.formal.cardio.intervention.contrast[[1,6]],5))),
                                                                ", ",bquote(.(round(secondary.str.g.moderator.formal.cardio.intervention.contrast[[1,7]],5))),"]")


secondary.str.g.moderator.formal.cardio.intervention.plot=ggplot(data = secondary.data.str%>%
                                                        mutate(formal.cardio.intervention=ifelse(formal.cardio.intervention=="No","Only Resistance Training","Concurrent Training")),
                                                      aes(x=avg.rir,
                                                          y=yi))+THEME+theme(legend.position = "right",
                                                                             plot.caption = element_text(face = "italic"))+
  facet_wrap(~formal.cardio.intervention)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Concurrent Training x <10 Estimated RIR",
       caption = secondary.str.g.moderator.formal.cardio.intervention.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = secondary.str.g.moderator.formal.cardio.intervention.emms%>%
                mutate(formal.cardio.intervention=ifelse(formal.cardio.intervention=="No","Only Resistance Training","Concurrent Training")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = secondary.str.g.moderator.formal.cardio.intervention.emms%>%
                mutate(formal.cardio.intervention=ifelse(formal.cardio.intervention=="No","Only Resistance Training","Concurrent Training")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = secondary.str.g.moderator.formal.cardio.intervention.emms%>%
              mutate(formal.cardio.intervention=ifelse(formal.cardio.intervention=="No","Only Resistance Training","Concurrent Training")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### g progressive overload
secondary.str.g.moderator.progression=secondary.str.g.best.model.prep%>%
  update(~. + log1p(avg.rir)*progression)%>%
  robust.rma.mv(cluster = secondary.data.str$study)

secondary.str.g.moderator.progression.emms=secondary.str.g.moderator.progression%>%emmprep(at=list(avg.rir=0:10))%>%
  emmeans(~avg.rir+progression,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.str.g.moderator.progression)

secondary.str.g.moderator.progression.contrast=emmprep(secondary.str.g.moderator.progression,
                                                           at=list(avg.rir=secondary.str.mean.rir))%>%
  emmeans(~avg.rir+progression,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()

secondary.str.g.moderator.progression.contrast.label=paste0('Interaction Contrast (No Overload - Progressive Overload) = ',
                                                                bquote(.(round(secondary.str.g.moderator.progression.contrast[[1,3]],5))),
                                                                " [95% CI: ",bquote(.(round(secondary.str.g.moderator.progression.contrast[[1,6]],5))),
                                                                ", ",bquote(.(round(secondary.str.g.moderator.progression.contrast[[1,7]],5))),"]")


secondary.str.g.moderator.progression.plot=ggplot(data = secondary.data.str%>%
                                                        mutate(progression=ifelse(progression=="No","No Progressive Overload","Progressive Overload")),
                                                      aes(x=avg.rir,
                                                          y=yi))+THEME+theme(legend.position = "right",
                                                                             plot.caption = element_text(face = "italic"))+
  facet_wrap(~progression)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Progressive Overload x <10 Estimated RIR",
       caption = secondary.str.g.moderator.progression.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = secondary.str.g.moderator.progression.emms%>%
                mutate(progression=ifelse(progression=="No","No Progressive Overload","Progressive Overload")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = secondary.str.g.moderator.progression.emms%>%
                mutate(progression=ifelse(progression=="No","No Progressive Overload","Progressive Overload")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = secondary.str.g.moderator.progression.emms%>%
              mutate(progression=ifelse(progression=="No","No Progressive Overload","Progressive Overload")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### g failure definition
secondary.str.g.moderator.failure.definition=secondary.str.g.best.model.prep%>%
  update(~. + log1p(avg.rir)*failure.definition)%>%
  robust.rma.mv(cluster = secondary.data.str$study)

secondary.str.g.moderator.failure.definition.emms=secondary.str.g.moderator.failure.definition%>%emmprep(at=list(avg.rir=0:10))%>%
  emmeans(~avg.rir+failure.definition,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.str.g.moderator.failure.definition)

secondary.str.g.moderator.failure.definition.contrast=emmprep(secondary.str.g.moderator.failure.definition,
                                                           at=list(avg.rir=secondary.str.mean.rir))%>%
  emmeans(~avg.rir+failure.definition,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()

secondary.str.g.moderator.failure.definition.contrast.label=paste0('Interaction Contrast (No Definition - Definition Provided) = ',
                                                                bquote(.(round(secondary.str.g.moderator.failure.definition.contrast[[1,3]],5))),
                                                                " [95% CI: ",bquote(.(round(secondary.str.g.moderator.failure.definition.contrast[[1,6]],5))),
                                                                ", ",bquote(.(round(secondary.str.g.moderator.failure.definition.contrast[[1,7]],5))),"]")


secondary.str.g.moderator.failure.definition.plot=ggplot(data = secondary.data.str%>%
                                                        mutate(failure.definition=ifelse(failure.definition=="No","No Failure Definition","Failure Definition Provided")),
                                                      aes(x=avg.rir,
                                                          y=yi))+THEME+theme(legend.position = "right",
                                                                             plot.caption = element_text(face = "italic"))+
  facet_wrap(~failure.definition)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Failure Definition x <10 Estimated RIR",
       caption = secondary.str.g.moderator.failure.definition.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = secondary.str.g.moderator.failure.definition.emms%>%
                mutate(failure.definition=ifelse(failure.definition=="No","No Failure Definition","Failure Definition Provided")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = secondary.str.g.moderator.failure.definition.emms%>%
                mutate(failure.definition=ifelse(failure.definition=="No","No Failure Definition","Failure Definition Provided")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = secondary.str.g.moderator.failure.definition.emms%>%
              mutate(failure.definition=ifelse(failure.definition=="No","No Failure Definition","Failure Definition Provided")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### g set structure
secondary.str.g.moderator.alternative.set.structure=secondary.str.g.best.model.prep%>%
  update(~. + log1p(avg.rir)*alternative.set.structure)%>%
  robust.rma.mv(cluster = secondary.data.str$study)

secondary.str.g.moderator.alternative.set.structure.emms=secondary.str.g.moderator.alternative.set.structure%>%emmprep(at=list(avg.rir=0:10))%>%
  emmeans(~avg.rir+alternative.set.structure,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.str.g.moderator.alternative.set.structure)

secondary.str.g.moderator.alternative.set.structure.contrast=emmprep(secondary.str.g.moderator.alternative.set.structure,
                                                           at=list(avg.rir=secondary.str.mean.rir))%>%
  emmeans(~avg.rir+alternative.set.structure,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()

secondary.str.g.moderator.alternative.set.structure.contrast.label=paste0('Interaction Contrast (Traditional - Alternative) = ',
                                                                bquote(.(round(secondary.str.g.moderator.alternative.set.structure.contrast[[1,3]],5))),
                                                                " [95% CI: ",bquote(.(round(secondary.str.g.moderator.alternative.set.structure.contrast[[1,6]],5))),
                                                                ", ",bquote(.(round(secondary.str.g.moderator.alternative.set.structure.contrast[[1,7]],5))),"]")


secondary.str.g.moderator.alternative.set.structure.plot=ggplot(data = secondary.data.str%>%
                                                        mutate(alternative.set.structure=ifelse(alternative.set.structure=="No","Traditional Set Structure","Alternative Set Structure")),
                                                      aes(x=avg.rir,
                                                          y=yi))+THEME+theme(legend.position = "right",
                                                                             plot.caption = element_text(face = "italic"))+
  facet_wrap(~alternative.set.structure)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Set Structure x <10 Estimated RIR",
       caption = secondary.str.g.moderator.alternative.set.structure.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = secondary.str.g.moderator.alternative.set.structure.emms%>%
                mutate(alternative.set.structure=ifelse(alternative.set.structure=="No","Traditional Set Structure","Alternative Set Structure")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = secondary.str.g.moderator.alternative.set.structure.emms%>%
                mutate(alternative.set.structure=ifelse(alternative.set.structure=="No","Traditional Set Structure","Alternative Set Structure")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = secondary.str.g.moderator.alternative.set.structure.emms%>%
              mutate(alternative.set.structure=ifelse(alternative.set.structure=="No","Traditional Set Structure","Alternative Set Structure")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### g intended velocity
secondary.str.g.moderator.train.intent=secondary.str.g.best.model.prep%>%
  update(~. + log1p(avg.rir)*train.intent)%>%
  robust.rma.mv(cluster = secondary.data.str$study)

secondary.str.g.moderator.train.intent.emms=secondary.str.g.moderator.train.intent%>%emmprep(at=list(avg.rir=0:10))%>%
  emmeans(~avg.rir+train.intent,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.str.g.moderator.train.intent)

secondary.str.g.moderator.train.intent.contrast=emmprep(secondary.str.g.moderator.train.intent,
                                                           at=list(avg.rir=secondary.str.mean.rir))%>%
  emmeans(~avg.rir+train.intent,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()

secondary.str.g.moderator.train.intent.contrast.label=paste0('Interaction Contrast (Submaximal - Maximal) = ',
                                                                bquote(.(round(secondary.str.g.moderator.train.intent.contrast[[1,3]],5))),
                                                                " [95% CI: ",bquote(.(round(secondary.str.g.moderator.train.intent.contrast[[1,6]],5))),
                                                                ", ",bquote(.(round(secondary.str.g.moderator.train.intent.contrast[[1,7]],5))),"]")


secondary.str.g.moderator.train.intent.plot=ggplot(data = secondary.data.str%>%
                                                        mutate(train.intent=ifelse(train.intent=="No","Submaximal Intended Concentric Velocity","Maximal Intended Concentric Velocity")),
                                                      aes(x=avg.rir,
                                                          y=yi))+THEME+theme(legend.position = "right",
                                                                             plot.caption = element_text(face = "italic"))+
  facet_wrap(~train.intent)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Intended Concentric Velocity x <10 Estimated RIR",
       caption = secondary.str.g.moderator.train.intent.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = secondary.str.g.moderator.train.intent.emms%>%
                mutate(train.intent=ifelse(train.intent=="No","Submaximal Intended Concentric Velocity","Maximal Intended Concentric Velocity")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = secondary.str.g.moderator.train.intent.emms%>%
                mutate(train.intent=ifelse(train.intent=="No","Submaximal Intended Concentric Velocity","Maximal Intended Concentric Velocity")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = secondary.str.g.moderator.train.intent.emms%>%
              mutate(train.intent=ifelse(train.intent=="No","Submaximal Intended Concentric Velocity","Maximal Intended Concentric Velocity")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### g type of strength test
secondary.str.g.moderator.isometric.isokinetic.isotonic=secondary.str.g.best.model.prep%>%
  update(~. + log1p(avg.rir)*isometric.isokinetic.isotonic)%>%
  robust.rma.mv(cluster = secondary.data.str$study)

secondary.str.g.moderator.isometric.isokinetic.isotonic.emms=secondary.str.g.moderator.isometric.isokinetic.isotonic%>%emmprep(at=list(avg.rir=0:10))%>%
  emmeans(~avg.rir+isometric.isokinetic.isotonic,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.str.g.moderator.isometric.isokinetic.isotonic)

secondary.str.g.moderator.isometric.isokinetic.isotonic.contrast=emmprep(secondary.str.g.moderator.isometric.isokinetic.isotonic,
                                                        at=list(avg.rir=secondary.str.mean.rir))%>%
  emmeans(~avg.rir+isometric.isokinetic.isotonic,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()

secondary.str.g.moderator.isometric.isokinetic.isotonic.contrast.label=paste0('Interaction Contrast (Isokinetic - Isometric) = ',
                                                             bquote(.(round(secondary.str.g.moderator.isometric.isokinetic.isotonic.contrast[[1,3]],5))),
                                                             " [95% CI: ",bquote(.(round(secondary.str.g.moderator.isometric.isokinetic.isotonic.contrast[[1,6]],5))),
                                                             ", ",bquote(.(round(secondary.str.g.moderator.isometric.isokinetic.isotonic.contrast[[1,7]],5))),"]\n",
                                                             'Interaction Contrast (Isokinetic - Isotonic) = ',
                                                             bquote(.(round(secondary.str.g.moderator.isometric.isokinetic.isotonic.contrast[[2,3]],5))),
                                                             " [95% CI: ",bquote(.(round(secondary.str.g.moderator.isometric.isokinetic.isotonic.contrast[[2,6]],5))),
                                                             ", ",bquote(.(round(secondary.str.g.moderator.isometric.isokinetic.isotonic.contrast[[2,7]],5))),"]\n",
                                                             'Interaction Contrast (Isometric - Isotonic) = ',
                                                             bquote(.(round(secondary.str.g.moderator.isometric.isokinetic.isotonic.contrast[[3,3]],5))),
                                                             " [95% CI: ",bquote(.(round(secondary.str.g.moderator.isometric.isokinetic.isotonic.contrast[[3,6]],5))),
                                                             ", ",bquote(.(round(secondary.str.g.moderator.isometric.isokinetic.isotonic.contrast[[3,7]],5))),"]")


secondary.str.g.moderator.isometric.isokinetic.isotonic.plot=ggplot(data = secondary.data.str%>%
                                                     mutate(isometric.isokinetic.isotonic=ifelse(isometric.isokinetic.isotonic=="isokinetic","Isokinetic",
                                                                                  ifelse(isometric.isokinetic.isotonic=="isometric","Isometric","Isotonic"))),
                                                   aes(x=avg.rir,
                                                       y=yi))+THEME+theme(legend.position = "right",
                                                                          plot.caption = element_text(face = "italic"))+
  facet_wrap(~isometric.isokinetic.isotonic)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Type of Strength Test x <10 Estimated RIR",
       caption = secondary.str.g.moderator.isometric.isokinetic.isotonic.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = secondary.str.g.moderator.isometric.isokinetic.isotonic.emms%>%
                mutate(isometric.isokinetic.isotonic=ifelse(isometric.isokinetic.isotonic=="isokinetic","Isokinetic",
                                                            ifelse(isometric.isokinetic.isotonic=="isometric","Isometric","Isotonic"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = secondary.str.g.moderator.isometric.isokinetic.isotonic.emms%>%
                mutate(isometric.isokinetic.isotonic=ifelse(isometric.isokinetic.isotonic=="isokinetic","Isokinetic",
                                                            ifelse(isometric.isokinetic.isotonic=="isometric","Isometric","Isotonic"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = secondary.str.g.moderator.isometric.isokinetic.isotonic.emms%>%
              mutate(isometric.isokinetic.isotonic=ifelse(isometric.isokinetic.isotonic=="isokinetic","Isokinetic",
                                                          ifelse(isometric.isokinetic.isotonic=="isometric","Isometric","Isotonic"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### g 1RM or non-1RM
secondary.str.g.moderator.RM.max.submax=secondary.str.g.best.model.prep%>%
  update(~. + log1p(avg.rir)*RM.max.submax)%>%
  robust.rma.mv(cluster = secondary.data.str$study)

secondary.str.g.moderator.RM.max.submax.emms=secondary.str.g.moderator.RM.max.submax%>%emmprep(at=list(avg.rir=0:10))%>%
  emmeans(~avg.rir+RM.max.submax,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.str.g.moderator.RM.max.submax)

secondary.str.g.moderator.RM.max.submax.contrast=emmprep(secondary.str.g.moderator.RM.max.submax,
                                                           at=list(avg.rir=secondary.str.mean.rir))%>%
  emmeans(~avg.rir+RM.max.submax,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()

secondary.str.g.moderator.RM.max.submax.contrast.label=paste0('Interaction Contrast (1RM - Non 1RM) = ',
                                                                bquote(.(round(secondary.str.g.moderator.RM.max.submax.contrast[[1,3]],5))),
                                                                " [95% CI: ",bquote(.(round(secondary.str.g.moderator.RM.max.submax.contrast[[1,6]],5))),
                                                                ", ",bquote(.(round(secondary.str.g.moderator.RM.max.submax.contrast[[1,7]],5))),"]")


secondary.str.g.moderator.RM.max.submax.plot=ggplot(data = secondary.data.str%>%
                                                        mutate(RM.max.submax=ifelse(RM.max.submax=="1RM","1RM","Non-1RM"))%>%
                                                    filter(RM.max.submax!="NA"),
                                                      aes(x=avg.rir,
                                                          y=yi))+THEME+theme(legend.position = "right",
                                                                             plot.caption = element_text(face = "italic"))+
  facet_wrap(~RM.max.submax)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Type of Isotonic Strength Test x <10 Estimated RIR",
       caption = secondary.str.g.moderator.RM.max.submax.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = secondary.str.g.moderator.RM.max.submax.emms%>%
                mutate(RM.max.submax=ifelse(RM.max.submax=="1RM","1RM","Non-1RM")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = secondary.str.g.moderator.RM.max.submax.emms%>%
                mutate(RM.max.submax=ifelse(RM.max.submax=="1RM","1RM","Non-1RM")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = secondary.str.g.moderator.RM.max.submax.emms%>%
              mutate(RM.max.submax=ifelse(RM.max.submax=="1RM","1RM","Non-1RM")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### combine g plots
secondary.str.g.moderator.doc=(secondary.str.g.moderator.load.plot/plot_spacer()/secondary.str.g.moderator.volume.method.plot/plot_spacer()/
  secondary.str.g.moderator.train.status.plot/plot_spacer()/secondary.str.g.moderator.weeks.plot/plot_spacer()/
  secondary.str.g.moderator.adj.sets.plot/plot_spacer()/secondary.str.g.moderator.frequency.per.muscle.plot/plot_spacer()/
  secondary.str.g.moderator.age.plot/plot_spacer()/secondary.str.g.moderator.sex.percent.male.plot/plot_spacer()/
  secondary.str.g.moderator.within.between.design.plot/plot_spacer()/secondary.str.g.moderator.upper.lower.other.plot/plot_spacer()/
  secondary.str.g.moderator.train.exercise.plot/plot_spacer()/secondary.str.g.moderator.formal.cardio.intervention.plot/plot_spacer()/
  secondary.str.g.moderator.progression.plot/plot_spacer()/secondary.str.g.moderator.failure.definition.plot/plot_spacer()/
  secondary.str.g.moderator.alternative.set.structure.plot/plot_spacer()/secondary.str.g.moderator.train.intent.plot/plot_spacer()/
  secondary.str.g.moderator.isometric.isokinetic.isotonic.plot/plot_spacer()/secondary.str.g.moderator.RM.max.submax.plot)&
  plot_annotation(title = "Moderator Analyses from Linear-log Multilevel Model for Maximal Strength (Standardized Mean Change)")&
  THEME

ggsave(secondary.str.g.moderator.doc,device = "pdf",filename="secondary.str.g.moderator.doc.pdf",width = 12,height = 156,limitsize = FALSE,
       path = "/Users/zacrobinson/Documents/For Publication PTF Meta/Plots")


#### rom load
secondary.str.rom.moderator.load=secondary.str.rom.best.model.prep%>%
  update(~. + avg.rir:load.set)%>%
  robust.rma.mv(cluster = secondary.data.str$study)

secondary.str.rom.moderator.load.emms=secondary.str.rom.moderator.load%>%
  emmprep(at=list(avg.rir=0:10,
                  load.set=c(30,60,90)))%>%
  emmeans(~avg.rir+load.set,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.str.rom.moderator.load)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

secondary.str.rom.moderator.load.contrast=emmprep(secondary.str.rom.moderator.load,
                                              at=list(avg.rir=secondary.str.mean.rir,
                                                      load.set=secondary.str.mean.load.set))%>%
  emmeans(~avg.rir+load.set,
          weights = "prop")%>%
  contrast("consec",interaction=TRUE)%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

secondary.str.rom.moderator.load.contrast.label=paste0('Interaction Contrast (Load) = ',
                                                   bquote(.(round(secondary.str.rom.moderator.load.contrast[[3]],5))),
                                                   " [95% CI: ",bquote(.(round(secondary.str.rom.moderator.load.contrast[[6]],5))),
                                                   ", ",bquote(.(round(secondary.str.rom.moderator.load.contrast[[7]],5))),"]")

secondary.str.rom.moderator.load.plot=ggplot(data = secondary.data.str%>%
                                           mutate(load.set=ifelse(load.set<45,"30% of 1RM",
                                                                  ifelse(load.set>45&load.set<75,"60% of 1RM","90% of 1RM"))),
                                         aes(x=avg.rir,
                                             y=yi.rom.exp))+THEME+theme(plot.caption = element_text(face = "italic"))+
  facet_wrap(~load.set)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Load (% of 1RM) x <10 Estimated RIR",
       caption = secondary.str.rom.moderator.load.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = secondary.str.rom.moderator.load.emms%>%
                mutate(load.set=ifelse(load.set==30,"30% of 1RM",
                                       ifelse(load.set==60,"60% of 1RM","90% of 1RM"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = secondary.str.rom.moderator.load.emms%>%
                mutate(load.set=ifelse(load.set==30,"30% of 1RM",
                                       ifelse(load.set==60,"60% of 1RM","90% of 1RM"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = secondary.str.rom.moderator.load.emms%>%
              mutate(load.set=ifelse(load.set==30,"30% of 1RM",
                                     ifelse(load.set==60,"60% of 1RM","90% of 1RM"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### rom volume equating method
secondary.str.rom.moderator.volume.method=secondary.str.rom.best.model.prep%>%
  update(~. + avg.rir:set.rep.equated)%>%
  robust.rma.mv(cluster = secondary.data.str$study)

secondary.str.rom.moderator.volume.method.emms=secondary.str.rom.moderator.volume.method%>%emmprep(at=list(avg.rir=0:10))%>%
  emmeans(~avg.rir+set.rep.equated,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.str.rom.moderator.volume.method)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

secondary.str.rom.moderator.volume.method.contrast=emmprep(secondary.str.rom.moderator.volume.method,
                                                       at=list(avg.rir=secondary.str.mean.rir))%>%
  emmeans(~avg.rir+set.rep.equated,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

secondary.str.rom.moderator.volume.method.contrast.label=paste0('Interaction Contrast (Set - Rep) = ',
                                                            bquote(.(round(secondary.str.rom.moderator.volume.method.contrast[[1,3]],5))),
                                                            " [95% CI: ",bquote(.(round(secondary.str.rom.moderator.volume.method.contrast[[1,6]],5))),
                                                            ", ",bquote(.(round(secondary.str.rom.moderator.volume.method.contrast[[1,7]],5))),"]\n",
                                                            'Interaction Contrast (Set - Both) = ',
                                                            bquote(.(round(secondary.str.rom.moderator.volume.method.contrast[[2,3]],5))),
                                                            " [95% CI: ",bquote(.(round(secondary.str.rom.moderator.volume.method.contrast[[2,6]],5))),
                                                            ", ",bquote(.(round(secondary.str.rom.moderator.volume.method.contrast[[2,7]],5))),"]\n",
                                                            'Interaction Contrast (Rep - Both) = ',
                                                            bquote(.(round(secondary.str.rom.moderator.volume.method.contrast[[3,3]],5))),
                                                            " [95% CI: ",bquote(.(round(secondary.str.rom.moderator.volume.method.contrast[[3,6]],5))),
                                                            ", ",bquote(.(round(secondary.str.rom.moderator.volume.method.contrast[[3,7]],5))),"]")


secondary.str.rom.moderator.volume.method.plot=ggplot(data = secondary.data.str%>%
                                                    mutate(set.rep.equated=ifelse(set.rep.equated=="set","Set Equated",
                                                                                  ifelse(set.rep.equated=="rep","Repetition Equated","Both Equated"))),
                                                  aes(x=avg.rir,
                                                      y=yi.rom.exp))+THEME+theme(legend.position = "right",
                                                                         plot.caption = element_text(face = "italic"))+
  facet_wrap(~set.rep.equated)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Method of Volume Equating x <10 Estimated RIR",
       caption = secondary.str.rom.moderator.volume.method.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = secondary.str.rom.moderator.volume.method.emms%>%
                mutate(set.rep.equated=ifelse(set.rep.equated=="set","Set Equated",
                                              ifelse(set.rep.equated=="rep","Repetition Equated","Both Equated"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = secondary.str.rom.moderator.volume.method.emms%>%
                mutate(set.rep.equated=ifelse(set.rep.equated=="set","Set Equated",
                                              ifelse(set.rep.equated=="rep","Repetition Equated","Both Equated"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = secondary.str.rom.moderator.volume.method.emms%>%
              mutate(set.rep.equated=ifelse(set.rep.equated=="set","Set Equated",
                                            ifelse(set.rep.equated=="rep","Repetition Equated","Both Equated"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### rom training status
secondary.str.rom.moderator.train.status=secondary.str.rom.best.model.prep%>%
  update(~. + avg.rir:train.status)%>%
  robust.rma.mv(cluster = secondary.data.str$study)

secondary.str.rom.moderator.train.status.emms=secondary.str.rom.moderator.train.status%>%emmprep(at=list(avg.rir=0:10))%>%
  emmeans(~avg.rir+train.status,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.str.rom.moderator.train.status)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

secondary.str.rom.moderator.train.status.contrast=emmprep(secondary.str.rom.moderator.train.status,
                                                      at=list(avg.rir=secondary.str.mean.rir))%>%
  emmeans(~avg.rir+train.status,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

secondary.str.rom.moderator.train.status.contrast.label=paste0('Interaction Contrast (Trained - Untrained) = ',
                                                           bquote(.(round(secondary.str.rom.moderator.train.status.contrast[[1,3]],5))),
                                                           " [95% CI: ",bquote(.(round(secondary.str.rom.moderator.train.status.contrast[[1,6]],5))),
                                                           ", ",bquote(.(round(secondary.str.rom.moderator.train.status.contrast[[1,7]],5))),"]")


secondary.str.rom.moderator.train.status.plot=ggplot(data = secondary.data.str%>%
                                                   mutate(train.status=ifelse(train.status=="trained","Trained","Untrained")),
                                                 aes(x=avg.rir,
                                                     y=yi.rom.exp))+THEME+theme(legend.position = "right",
                                                                        plot.caption = element_text(face = "italic"))+
  facet_wrap(~train.status)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Training Status x <10 Estimated RIR",
       caption = secondary.str.rom.moderator.train.status.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = secondary.str.rom.moderator.train.status.emms%>%
                mutate(train.status=ifelse(train.status=="trained","Trained","Untrained")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = secondary.str.rom.moderator.train.status.emms%>%
                mutate(train.status=ifelse(train.status=="trained","Trained","Untrained")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = secondary.str.rom.moderator.train.status.emms%>%
              mutate(train.status=ifelse(train.status=="trained","Trained","Untrained")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### rom weeks
secondary.str.rom.moderator.weeks=secondary.str.rom.best.model.prep%>%
  update(~. + avg.rir:weeks)%>%
  robust.rma.mv(cluster = secondary.data.str$study)

secondary.str.rom.moderator.weeks.emms=secondary.str.rom.moderator.weeks%>%
  emmprep(at=list(avg.rir=0:10,
                  weeks=c(4,8,12)))%>%
  emmeans(~avg.rir+weeks,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.str.rom.moderator.weeks)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

secondary.str.rom.moderator.weeks.contrast=emmprep(secondary.str.rom.moderator.weeks,
                                               at=list(avg.rir=secondary.str.mean.rir,
                                                       weeks=secondary.str.mean.weeks))%>%
  emmeans(~avg.rir+weeks,
          weights = "prop")%>%
  contrast("consec",interaction=TRUE)%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

secondary.str.rom.moderator.weeks.contrast.label=paste0('Interaction Contrast (Weeks) = ',
                                                    bquote(.(round(secondary.str.rom.moderator.weeks.contrast[[3]],5))),
                                                    " [95% CI: ",bquote(.(round(secondary.str.rom.moderator.weeks.contrast[[6]],5))),
                                                    ", ",bquote(.(round(secondary.str.rom.moderator.weeks.contrast[[7]],5))),"]")

secondary.str.rom.moderator.weeks.plot=ggplot(data = secondary.data.str%>%
                                            mutate(weeks=ifelse(weeks<6,"4 Weeks",
                                                                ifelse(weeks>6&weeks<10,"8 Weeks","12 Weeks"))),
                                          aes(x=avg.rir,
                                              y=yi.rom.exp))+THEME+theme(plot.caption = element_text(face = "italic"))+
  facet_wrap(~fct_relevel(weeks,c("4 Weeks","8 Weeks","12 Weeks")))+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Intervention Duration (Weeks) x <10 Estimated RIR",
       caption = secondary.str.rom.moderator.weeks.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = secondary.str.rom.moderator.weeks.emms%>%
                mutate(weeks=ifelse(weeks==4,"4 Weeks",
                                    ifelse(weeks==8,"8 Weeks","12 Weeks"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = secondary.str.rom.moderator.weeks.emms%>%
                mutate(weeks=ifelse(weeks==4,"4 Weeks",
                                    ifelse(weeks==8,"8 Weeks","12 Weeks"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = secondary.str.rom.moderator.weeks.emms%>%
              mutate(weeks=ifelse(weeks==4,"4 Weeks",
                                  ifelse(weeks==8,"8 Weeks","12 Weeks"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### rom adjusted sets 
secondary.str.rom.moderator.adj.sets=secondary.str.rom.best.model.prep%>%
  update(~. + avg.rir*adj.sets.week)%>%
  robust.rma.mv(cluster = secondary.data.str$study)

secondary.str.rom.moderator.adj.sets.emms=secondary.str.rom.moderator.adj.sets%>%
  emmprep(at=list(avg.rir=0:10,
                  adj.sets.week=c(8,16,24)))%>%
  emmeans(~avg.rir+adj.sets.week,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.str.rom.moderator.adj.sets)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

secondary.str.rom.moderator.adj.sets.contrast=emmprep(secondary.str.rom.moderator.adj.sets,
                                                  at=list(avg.rir=secondary.str.mean.rir,
                                                          adj.sets.week=secondary.str.mean.adj.sets.week))%>%
  emmeans(~avg.rir+adj.sets.week,
          weights = "prop")%>%
  contrast("consec",interaction=TRUE)%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

secondary.str.rom.moderator.adj.sets.contrast.label=paste0('Interaction Contrast (Adjusted Direct Sets) = ',
                                                       bquote(.(round(secondary.str.rom.moderator.adj.sets.contrast[[3]],5))),
                                                       " [95% CI: ",bquote(.(round(secondary.str.rom.moderator.adj.sets.contrast[[6]],5))),
                                                       ", ",bquote(.(round(secondary.str.rom.moderator.adj.sets.contrast[[7]],5))),"]")

secondary.str.rom.moderator.adj.sets.plot=ggplot(data = secondary.data.str%>%
                                               mutate(adj.sets.week=ifelse(adj.sets.week<12,"8 Sets",
                                                                           ifelse(adj.sets.week>12&adj.sets.week<20,"16 Sets","24 Sets"))),
                                             aes(x=avg.rir,
                                                 y=yi.rom.exp))+THEME+theme(plot.caption = element_text(face = "italic"))+
  facet_wrap(~fct_relevel(adj.sets.week,c("8 Sets","16 Sets","24 Sets")))+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Adjusted Direct Sets Per Week Per Muscle Group x <10 Estimated RIR",
       caption = secondary.str.rom.moderator.adj.sets.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = secondary.str.rom.moderator.adj.sets.emms%>%
                mutate(adj.sets.week=ifelse(adj.sets.week==8,"8 Sets",
                                            ifelse(adj.sets.week==16,"16 Sets","24 Sets"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = secondary.str.rom.moderator.adj.sets.emms%>%
                mutate(adj.sets.week=ifelse(adj.sets.week==8,"8 Sets",
                                            ifelse(adj.sets.week==16,"16 Sets","24 Sets"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = secondary.str.rom.moderator.adj.sets.emms%>%
              mutate(adj.sets.week=ifelse(adj.sets.week==8,"8 Sets",
                                          ifelse(adj.sets.week==16,"16 Sets","24 Sets"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### rom frequency
secondary.str.rom.moderator.frequency.per.muscle=secondary.str.rom.best.model.prep%>%
  update(~. + avg.rir*frequency.per.muscle)%>%
  robust.rma.mv(cluster = secondary.data.str$study)

secondary.str.rom.moderator.frequency.per.muscle.emms=secondary.str.rom.moderator.frequency.per.muscle%>%
  emmprep(at=list(avg.rir=0:10,
                  frequency.per.muscle=c(1,2,3)))%>%
  emmeans(~avg.rir+frequency.per.muscle,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.str.rom.moderator.frequency.per.muscle)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

secondary.str.rom.moderator.frequency.per.muscle.contrast=emmprep(secondary.str.rom.moderator.frequency.per.muscle,
                                                              at=list(avg.rir=secondary.str.mean.rir,
                                                                      frequency.per.muscle=secondary.str.mean.frequency.per.muscle))%>%
  emmeans(~avg.rir+frequency.per.muscle,
          weights = "prop")%>%
  contrast("consec",interaction=TRUE)%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

secondary.str.rom.moderator.frequency.per.muscle.contrast.label=paste0('Interaction Contrast (Direct Frequency) = ',
                                                                   bquote(.(round(secondary.str.rom.moderator.frequency.per.muscle.contrast[[3]],5))),
                                                                   " [95% CI: ",bquote(.(round(secondary.str.rom.moderator.frequency.per.muscle.contrast[[6]],5))),
                                                                   ", ",bquote(.(round(secondary.str.rom.moderator.frequency.per.muscle.contrast[[7]],5))),"]")

secondary.str.rom.moderator.frequency.per.muscle.plot=ggplot(data = secondary.data.str%>%
                                                           mutate(frequency.per.muscle=ifelse(frequency.per.muscle<1.5,"1x Per Week",
                                                                                              ifelse(frequency.per.muscle>1.5&frequency.per.muscle<2.5,"2x Per Week","3x Per Week"))),
                                                         aes(x=avg.rir,
                                                             y=yi.rom.exp))+THEME+theme(plot.caption = element_text(face = "italic"))+
  facet_wrap(~frequency.per.muscle)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Direct Frequency Per Week Per Muscle Group x <10 Estimated RIR",
       caption = secondary.str.rom.moderator.frequency.per.muscle.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = secondary.str.rom.moderator.frequency.per.muscle.emms%>%
                mutate(frequency.per.muscle=ifelse(frequency.per.muscle==1,"1x Per Week",
                                                   ifelse(frequency.per.muscle==2,"2x Per Week","3x Per Week"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = secondary.str.rom.moderator.frequency.per.muscle.emms%>%
                mutate(frequency.per.muscle=ifelse(frequency.per.muscle==1,"1x Per Week",
                                                   ifelse(frequency.per.muscle==2,"2x Per Week","3x Per Week"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = secondary.str.rom.moderator.frequency.per.muscle.emms%>%
              mutate(frequency.per.muscle=ifelse(frequency.per.muscle==1,"1x Per Week",
                                                 ifelse(frequency.per.muscle==2,"2x Per Week","3x Per Week"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### rom age
secondary.str.rom.moderator.age=secondary.str.rom.best.model.prep%>%
  update(~. + avg.rir*age)%>%
  robust.rma.mv(cluster = secondary.data.str$study)

secondary.str.rom.moderator.age.emms=secondary.str.rom.moderator.age%>%
  emmprep(at=list(avg.rir=0:10,
                  age=c(20,40,60)))%>%
  emmeans(~avg.rir+age,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.str.rom.moderator.age)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

secondary.str.rom.moderator.age.contrast=emmprep(secondary.str.rom.moderator.age,
                                             at=list(avg.rir=secondary.str.mean.rir,
                                                     age=secondary.str.mean.age))%>%
  emmeans(~avg.rir+age,
          weights = "prop")%>%
  contrast("consec",interaction=TRUE)%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

secondary.str.rom.moderator.age.contrast.label=paste0('Interaction Contrast (Age) = ',
                                                  bquote(.(round(secondary.str.rom.moderator.age.contrast[[3]],5))),
                                                  " [95% CI: ",bquote(.(round(secondary.str.rom.moderator.age.contrast[[6]],5))),
                                                  ", ",bquote(.(round(secondary.str.rom.moderator.age.contrast[[7]],5))),"]")

secondary.str.rom.moderator.age.plot=ggplot(data = secondary.data.str%>%
                                          mutate(age=ifelse(age<30,"20 Years Old",
                                                            ifelse(age>30&age<50,"40 Years Old","60 Years Old"))),
                                        aes(x=avg.rir,
                                            y=yi.rom.exp))+THEME+theme(plot.caption = element_text(face = "italic"))+
  facet_wrap(~age)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Age x <10 Estimated RIR",
       caption = secondary.str.rom.moderator.age.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = secondary.str.rom.moderator.age.emms%>%
                mutate(age=ifelse(age==20,"20 Years Old",
                                  ifelse(age==40,"40 Years Old","60 Years Old"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = secondary.str.rom.moderator.age.emms%>%
                mutate(age=ifelse(age==20,"20 Years Old",
                                  ifelse(age==40,"40 Years Old","60 Years Old"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = secondary.str.rom.moderator.age.emms%>%
              mutate(age=ifelse(age==20,"20 Years Old",
                                ifelse(age==40,"40 Years Old","60 Years Old"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### rom sex
secondary.str.rom.moderator.sex.percent.male=secondary.str.rom.best.model.prep%>%
  update(~. + avg.rir*sex.percent.male)%>%
  robust.rma.mv(cluster = secondary.data.str$study)

secondary.str.rom.moderator.sex.percent.male.emms=secondary.str.rom.moderator.sex.percent.male%>%
  emmprep(at=list(avg.rir=0:10,
                  sex.percent.male=c(0,100)))%>%
  emmeans(~avg.rir+sex.percent.male,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.str.rom.moderator.sex.percent.male)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

secondary.str.rom.moderator.sex.percent.male.contrast=emmprep(secondary.str.rom.moderator.sex.percent.male,
                                                          at=list(avg.rir=secondary.str.mean.rir,
                                                                  sex.percent.male=secondary.str.mean.sex))%>%
  emmeans(~avg.rir+sex.percent.male,
          weights = "prop")%>%
  contrast("consec",interaction=TRUE)%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

secondary.str.rom.moderator.sex.percent.male.contrast.label=paste0('Interaction Contrast (Sex) = ',
                                                               bquote(.(round(secondary.str.rom.moderator.sex.percent.male.contrast[[3]],5))),
                                                               " [95% CI: ",bquote(.(round(secondary.str.rom.moderator.sex.percent.male.contrast[[6]],5))),
                                                               ", ",bquote(.(round(secondary.str.rom.moderator.sex.percent.male.contrast[[7]],5))),"]")

secondary.str.rom.moderator.sex.percent.male.plot=ggplot(data = secondary.data.str%>%
                                                       mutate(sex.percent.male=ifelse(sex.percent.male<50,"Female","Male")),
                                                     aes(x=avg.rir,
                                                         y=yi.rom.exp))+THEME+theme(plot.caption = element_text(face = "italic"))+
  facet_wrap(~sex.percent.male)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Sex (% Male) x <10 Estimated RIR",
       caption = secondary.str.rom.moderator.sex.percent.male.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = secondary.str.rom.moderator.sex.percent.male.emms%>%
                mutate(sex.percent.male=ifelse(sex.percent.male==0,"Female","Male")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = secondary.str.rom.moderator.sex.percent.male.emms%>%
                mutate(sex.percent.male=ifelse(sex.percent.male==0,"Female","Male")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = secondary.str.rom.moderator.sex.percent.male.emms%>%
              mutate(sex.percent.male=ifelse(sex.percent.male==0,"Female","Male")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### rom study design
secondary.str.rom.moderator.within.between.design=secondary.str.rom.best.model.prep%>%
  update(~. + avg.rir*within.between.design)%>%
  robust.rma.mv(cluster = secondary.data.str$study)

secondary.str.rom.moderator.within.between.design.emms=secondary.str.rom.moderator.within.between.design%>%emmprep(at=list(avg.rir=0:10))%>%
  emmeans(~avg.rir+within.between.design,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.str.rom.moderator.within.between.design)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

secondary.str.rom.moderator.within.between.design.contrast=emmprep(secondary.str.rom.moderator.within.between.design,
                                                               at=list(avg.rir=secondary.str.mean.rir))%>%
  emmeans(~avg.rir+within.between.design,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

secondary.str.rom.moderator.within.between.design.contrast.label=paste0('Interaction Contrast (Between - Within) = ',
                                                                    bquote(.(round(secondary.str.rom.moderator.within.between.design.contrast[[1,3]],5))),
                                                                    " [95% CI: ",bquote(.(round(secondary.str.rom.moderator.within.between.design.contrast[[1,6]],5))),
                                                                    ", ",bquote(.(round(secondary.str.rom.moderator.within.between.design.contrast[[1,7]],5))),"]")


secondary.str.rom.moderator.within.between.design.plot=ggplot(data = secondary.data.str%>%
                                                            mutate(within.between.design=ifelse(within.between.design=="Between","Between Subjects","Within Subjects")),
                                                          aes(x=avg.rir,
                                                              y=yi.rom.exp))+THEME+theme(legend.position = "right",
                                                                                 plot.caption = element_text(face = "italic"))+
  facet_wrap(~within.between.design)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Study Design x <10 Estimated RIR",
       caption = secondary.str.rom.moderator.within.between.design.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = secondary.str.rom.moderator.within.between.design.emms%>%
                mutate(within.between.design=ifelse(within.between.design=="Between","Between Subjects","Within Subjects")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = secondary.str.rom.moderator.within.between.design.emms%>%
                mutate(within.between.design=ifelse(within.between.design=="Between","Between Subjects","Within Subjects")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = secondary.str.rom.moderator.within.between.design.emms%>%
              mutate(within.between.design=ifelse(within.between.design=="Between","Between Subjects","Within Subjects")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### rom upper / lower outcome
secondary.str.rom.moderator.upper.lower.other=secondary.str.rom.best.model.prep%>%
  update(~. + avg.rir*upper.lower.other)%>%
  robust.rma.mv(cluster = secondary.data.str$study)

secondary.str.rom.moderator.upper.lower.other.emms=secondary.str.rom.moderator.upper.lower.other%>%emmprep(at=list(avg.rir=0:10))%>%
  emmeans(~avg.rir+upper.lower.other,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.str.rom.moderator.upper.lower.other)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

secondary.str.rom.moderator.upper.lower.other.contrast=emmprep(secondary.str.rom.moderator.upper.lower.other,
                                                           at=list(avg.rir=secondary.str.mean.rir))%>%
  emmeans(~avg.rir+upper.lower.other,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

secondary.str.rom.moderator.upper.lower.other.contrast.label=paste0('Interaction Contrast (Lower - Upper) = ',
                                                                bquote(.(round(secondary.str.rom.moderator.upper.lower.other.contrast[[1,3]],5))),
                                                                " [95% CI: ",bquote(.(round(secondary.str.rom.moderator.upper.lower.other.contrast[[1,6]],5))),
                                                                ", ",bquote(.(round(secondary.str.rom.moderator.upper.lower.other.contrast[[1,7]],5))),"]")


secondary.str.rom.moderator.upper.lower.other.plot=ggplot(data = secondary.data.str%>%
                                                        mutate(upper.lower.other=ifelse(upper.lower.other=="Lower","Lower Body Outcomes","Upper Body Outcomes")),
                                                      aes(x=avg.rir,
                                                          y=yi.rom.exp))+THEME+theme(legend.position = "right",
                                                                             plot.caption = element_text(face = "italic"))+
  facet_wrap(~upper.lower.other)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Upper/Lower Body Outcome x <10 Estimated RIR",
       caption = secondary.str.rom.moderator.upper.lower.other.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = secondary.str.rom.moderator.upper.lower.other.emms%>%
                mutate(upper.lower.other=ifelse(upper.lower.other=="Lower","Lower Body Outcomes","Upper Body Outcomes")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = secondary.str.rom.moderator.upper.lower.other.emms%>%
                mutate(upper.lower.other=ifelse(upper.lower.other=="Lower","Lower Body Outcomes","Upper Body Outcomes")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = secondary.str.rom.moderator.upper.lower.other.emms%>%
              mutate(upper.lower.other=ifelse(upper.lower.other=="Lower","Lower Body Outcomes","Upper Body Outcomes")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### rom exercise selection
secondary.str.rom.moderator.train.exercise=secondary.str.rom.best.model.prep%>%
  update(~. + avg.rir*train.exercise)%>%
  robust.rma.mv(cluster = secondary.data.str$study)

secondary.str.rom.moderator.train.exercise.emms=secondary.str.rom.moderator.train.exercise%>%emmprep(at=list(avg.rir=0:10))%>%
  emmeans(~avg.rir+train.exercise,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.str.rom.moderator.train.exercise)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

secondary.str.rom.moderator.train.exercise.contrast=emmprep(secondary.str.rom.moderator.train.exercise,
                                                        at=list(avg.rir=secondary.str.mean.rir))%>%
  emmeans(~avg.rir+train.exercise,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

secondary.str.rom.moderator.train.exercise.contrast.label=paste0('Interaction Contrast (Both - Multi-joint) = ',
                                                             bquote(.(round(secondary.str.rom.moderator.train.exercise.contrast[[1,3]],5))),
                                                             " [95% CI: ",bquote(.(round(secondary.str.rom.moderator.train.exercise.contrast[[1,6]],5))),
                                                             ", ",bquote(.(round(secondary.str.rom.moderator.train.exercise.contrast[[1,7]],5))),"]\n",
                                                             'Interaction Contrast (Both - Single-joint) = ',
                                                             bquote(.(round(secondary.str.rom.moderator.train.exercise.contrast[[2,3]],5))),
                                                             " [95% CI: ",bquote(.(round(secondary.str.rom.moderator.train.exercise.contrast[[2,6]],5))),
                                                             ", ",bquote(.(round(secondary.str.rom.moderator.train.exercise.contrast[[2,7]],5))),"]\n",
                                                             'Interaction Contrast (Multi-joint - Single-joint) = ',
                                                             bquote(.(round(secondary.str.rom.moderator.train.exercise.contrast[[3,3]],5))),
                                                             " [95% CI: ",bquote(.(round(secondary.str.rom.moderator.train.exercise.contrast[[3,6]],5))),
                                                             ", ",bquote(.(round(secondary.str.rom.moderator.train.exercise.contrast[[3,7]],5))),"]")


secondary.str.rom.moderator.train.exercise.plot=ggplot(data = secondary.data.str%>%
                                                     mutate(train.exercise=ifelse(train.exercise=="Both","Both",
                                                                                  ifelse(train.exercise=="Multi-joint","Multi-joint Exercises","Single-joint Exercises"))),
                                                   aes(x=avg.rir,
                                                       y=yi.rom.exp))+THEME+theme(legend.position = "right",
                                                                          plot.caption = element_text(face = "italic"))+
  facet_wrap(~train.exercise)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Exercise Selection x <10 Estimated RIR",
       caption = secondary.str.rom.moderator.train.exercise.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = secondary.str.rom.moderator.train.exercise.emms%>%
                mutate(train.exercise=ifelse(train.exercise=="Both","Both",
                                             ifelse(train.exercise=="Multi-joint","Multi-joint Exercises","Single-joint Exercises"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = secondary.str.rom.moderator.train.exercise.emms%>%
                mutate(train.exercise=ifelse(train.exercise=="Both","Both",
                                             ifelse(train.exercise=="Multi-joint","Multi-joint Exercises","Single-joint Exercises"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = secondary.str.rom.moderator.train.exercise.emms%>%
              mutate(train.exercise=ifelse(train.exercise=="Both","Both",
                                           ifelse(train.exercise=="Multi-joint","Multi-joint Exercises","Single-joint Exercises"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### rom concurrent training
secondary.str.rom.moderator.formal.cardio.intervention=secondary.str.rom.best.model.prep%>%
  update(~. + avg.rir*formal.cardio.intervention)%>%
  robust.rma.mv(cluster = secondary.data.str$study)

secondary.str.rom.moderator.formal.cardio.intervention.emms=secondary.str.rom.moderator.formal.cardio.intervention%>%emmprep(at=list(avg.rir=0:10))%>%
  emmeans(~avg.rir+formal.cardio.intervention,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.str.rom.moderator.formal.cardio.intervention)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

secondary.str.rom.moderator.formal.cardio.intervention.contrast=emmprep(secondary.str.rom.moderator.formal.cardio.intervention,
                                                                    at=list(avg.rir=secondary.str.mean.rir))%>%
  emmeans(~avg.rir+formal.cardio.intervention,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

secondary.str.rom.moderator.formal.cardio.intervention.contrast.label=paste0('Interaction Contrast (RT Only - Concurrent) = ',
                                                                         bquote(.(round(secondary.str.rom.moderator.formal.cardio.intervention.contrast[[1,3]],5))),
                                                                         " [95% CI: ",bquote(.(round(secondary.str.rom.moderator.formal.cardio.intervention.contrast[[1,6]],5))),
                                                                         ", ",bquote(.(round(secondary.str.rom.moderator.formal.cardio.intervention.contrast[[1,7]],5))),"]")


secondary.str.rom.moderator.formal.cardio.intervention.plot=ggplot(data = secondary.data.str%>%
                                                                 mutate(formal.cardio.intervention=ifelse(formal.cardio.intervention=="No","Only Resistance Training","Concurrent Training")),
                                                               aes(x=avg.rir,
                                                                   y=yi.rom.exp))+THEME+theme(legend.position = "right",
                                                                                      plot.caption = element_text(face = "italic"))+
  facet_wrap(~formal.cardio.intervention)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Concurrent Training x <10 Estimated RIR",
       caption = secondary.str.rom.moderator.formal.cardio.intervention.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = secondary.str.rom.moderator.formal.cardio.intervention.emms%>%
                mutate(formal.cardio.intervention=ifelse(formal.cardio.intervention=="No","Only Resistance Training","Concurrent Training")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = secondary.str.rom.moderator.formal.cardio.intervention.emms%>%
                mutate(formal.cardio.intervention=ifelse(formal.cardio.intervention=="No","Only Resistance Training","Concurrent Training")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = secondary.str.rom.moderator.formal.cardio.intervention.emms%>%
              mutate(formal.cardio.intervention=ifelse(formal.cardio.intervention=="No","Only Resistance Training","Concurrent Training")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### rom progressive overload
secondary.str.rom.moderator.progression=secondary.str.rom.best.model.prep%>%
  update(~. + avg.rir*progression)%>%
  robust.rma.mv(cluster = secondary.data.str$study)

secondary.str.rom.moderator.progression.emms=secondary.str.rom.moderator.progression%>%emmprep(at=list(avg.rir=0:10))%>%
  emmeans(~avg.rir+progression,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.str.rom.moderator.progression)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

secondary.str.rom.moderator.progression.contrast=emmprep(secondary.str.rom.moderator.progression,
                                                     at=list(avg.rir=secondary.str.mean.rir))%>%
  emmeans(~avg.rir+progression,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

secondary.str.rom.moderator.progression.contrast.label=paste0('Interaction Contrast (No Overload - Progressive Overload) = ',
                                                          bquote(.(round(secondary.str.rom.moderator.progression.contrast[[1,3]],5))),
                                                          " [95% CI: ",bquote(.(round(secondary.str.rom.moderator.progression.contrast[[1,6]],5))),
                                                          ", ",bquote(.(round(secondary.str.rom.moderator.progression.contrast[[1,7]],5))),"]")


secondary.str.rom.moderator.progression.plot=ggplot(data = secondary.data.str%>%
                                                  mutate(progression=ifelse(progression=="No","No Progressive Overload","Progressive Overload")),
                                                aes(x=avg.rir,
                                                    y=yi.rom.exp))+THEME+theme(legend.position = "right",
                                                                       plot.caption = element_text(face = "italic"))+
  facet_wrap(~progression)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Progressive Overload x <10 Estimated RIR",
       caption = secondary.str.rom.moderator.progression.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = secondary.str.rom.moderator.progression.emms%>%
                mutate(progression=ifelse(progression=="No","No Progressive Overload","Progressive Overload")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = secondary.str.rom.moderator.progression.emms%>%
                mutate(progression=ifelse(progression=="No","No Progressive Overload","Progressive Overload")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = secondary.str.rom.moderator.progression.emms%>%
              mutate(progression=ifelse(progression=="No","No Progressive Overload","Progressive Overload")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### rom failure definition
secondary.str.rom.moderator.failure.definition=secondary.str.rom.best.model.prep%>%
  update(~. + avg.rir*failure.definition)%>%
  robust.rma.mv(cluster = secondary.data.str$study)

secondary.str.rom.moderator.failure.definition.emms=secondary.str.rom.moderator.failure.definition%>%emmprep(at=list(avg.rir=0:10))%>%
  emmeans(~avg.rir+failure.definition,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.str.rom.moderator.failure.definition)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

secondary.str.rom.moderator.failure.definition.contrast=emmprep(secondary.str.rom.moderator.failure.definition,
                                                            at=list(avg.rir=secondary.str.mean.rir))%>%
  emmeans(~avg.rir+failure.definition,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

secondary.str.rom.moderator.failure.definition.contrast.label=paste0('Interaction Contrast (No Definition - Definition Provided) = ',
                                                                 bquote(.(round(secondary.str.rom.moderator.failure.definition.contrast[[1,3]],5))),
                                                                 " [95% CI: ",bquote(.(round(secondary.str.rom.moderator.failure.definition.contrast[[1,6]],5))),
                                                                 ", ",bquote(.(round(secondary.str.rom.moderator.failure.definition.contrast[[1,7]],5))),"]")


secondary.str.rom.moderator.failure.definition.plot=ggplot(data = secondary.data.str%>%
                                                         mutate(failure.definition=ifelse(failure.definition=="No","No Failure Definition","Failure Definition Provided")),
                                                       aes(x=avg.rir,
                                                           y=yi.rom.exp))+THEME+theme(legend.position = "right",
                                                                              plot.caption = element_text(face = "italic"))+
  facet_wrap(~failure.definition)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Failure Definition x <10 Estimated RIR",
       caption = secondary.str.rom.moderator.failure.definition.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = secondary.str.rom.moderator.failure.definition.emms%>%
                mutate(failure.definition=ifelse(failure.definition=="No","No Failure Definition","Failure Definition Provided")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = secondary.str.rom.moderator.failure.definition.emms%>%
                mutate(failure.definition=ifelse(failure.definition=="No","No Failure Definition","Failure Definition Provided")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = secondary.str.rom.moderator.failure.definition.emms%>%
              mutate(failure.definition=ifelse(failure.definition=="No","No Failure Definition","Failure Definition Provided")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### rom set structure
secondary.str.rom.moderator.alternative.set.structure=secondary.str.rom.best.model.prep%>%
  update(~. + avg.rir*alternative.set.structure)%>%
  robust.rma.mv(cluster = secondary.data.str$study)

secondary.str.rom.moderator.alternative.set.structure.emms=secondary.str.rom.moderator.alternative.set.structure%>%emmprep(at=list(avg.rir=0:10))%>%
  emmeans(~avg.rir+alternative.set.structure,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.str.rom.moderator.alternative.set.structure)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

secondary.str.rom.moderator.alternative.set.structure.contrast=emmprep(secondary.str.rom.moderator.alternative.set.structure,
                                                                   at=list(avg.rir=secondary.str.mean.rir))%>%
  emmeans(~avg.rir+alternative.set.structure,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

secondary.str.rom.moderator.alternative.set.structure.contrast.label=paste0('Interaction Contrast (Traditional - Alternative) = ',
                                                                        bquote(.(round(secondary.str.rom.moderator.alternative.set.structure.contrast[[1,3]],5))),
                                                                        " [95% CI: ",bquote(.(round(secondary.str.rom.moderator.alternative.set.structure.contrast[[1,6]],5))),
                                                                        ", ",bquote(.(round(secondary.str.rom.moderator.alternative.set.structure.contrast[[1,7]],5))),"]")


secondary.str.rom.moderator.alternative.set.structure.plot=ggplot(data = secondary.data.str%>%
                                                                mutate(alternative.set.structure=ifelse(alternative.set.structure=="No","Traditional Set Structure","Alternative Set Structure")),
                                                              aes(x=avg.rir,
                                                                  y=yi.rom.exp))+THEME+theme(legend.position = "right",
                                                                                     plot.caption = element_text(face = "italic"))+
  facet_wrap(~alternative.set.structure)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Set Structure x <10 Estimated RIR",
       caption = secondary.str.rom.moderator.alternative.set.structure.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = secondary.str.rom.moderator.alternative.set.structure.emms%>%
                mutate(alternative.set.structure=ifelse(alternative.set.structure=="No","Traditional Set Structure","Alternative Set Structure")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = secondary.str.rom.moderator.alternative.set.structure.emms%>%
                mutate(alternative.set.structure=ifelse(alternative.set.structure=="No","Traditional Set Structure","Alternative Set Structure")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = secondary.str.rom.moderator.alternative.set.structure.emms%>%
              mutate(alternative.set.structure=ifelse(alternative.set.structure=="No","Traditional Set Structure","Alternative Set Structure")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### rom intended velocity
secondary.str.rom.moderator.train.intent=secondary.str.rom.best.model.prep%>%
  update(~. + avg.rir*train.intent)%>%
  robust.rma.mv(cluster = secondary.data.str$study)

secondary.str.rom.moderator.train.intent.emms=secondary.str.rom.moderator.train.intent%>%emmprep(at=list(avg.rir=0:10))%>%
  emmeans(~avg.rir+train.intent,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.str.rom.moderator.train.intent)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

secondary.str.rom.moderator.train.intent.contrast=emmprep(secondary.str.rom.moderator.train.intent,
                                                      at=list(avg.rir=secondary.str.mean.rir))%>%
  emmeans(~avg.rir+train.intent,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

secondary.str.rom.moderator.train.intent.contrast.label=paste0('Interaction Contrast (Submaximal - Maximal) = ',
                                                           bquote(.(round(secondary.str.rom.moderator.train.intent.contrast[[1,3]],5))),
                                                           " [95% CI: ",bquote(.(round(secondary.str.rom.moderator.train.intent.contrast[[1,6]],5))),
                                                           ", ",bquote(.(round(secondary.str.rom.moderator.train.intent.contrast[[1,7]],5))),"]")


secondary.str.rom.moderator.train.intent.plot=ggplot(data = secondary.data.str%>%
                                                   mutate(train.intent=ifelse(train.intent=="No","Submaximal Intended Concentric Velocity","Maximal Intended Concentric Velocity")),
                                                 aes(x=avg.rir,
                                                     y=yi.rom.exp))+THEME+theme(legend.position = "right",
                                                                        plot.caption = element_text(face = "italic"))+
  facet_wrap(~train.intent)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Intended Concentric Velocity x <10 Estimated RIR",
       caption = secondary.str.rom.moderator.train.intent.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = secondary.str.rom.moderator.train.intent.emms%>%
                mutate(train.intent=ifelse(train.intent=="No","Submaximal Intended Concentric Velocity","Maximal Intended Concentric Velocity")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = secondary.str.rom.moderator.train.intent.emms%>%
                mutate(train.intent=ifelse(train.intent=="No","Submaximal Intended Concentric Velocity","Maximal Intended Concentric Velocity")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = secondary.str.rom.moderator.train.intent.emms%>%
              mutate(train.intent=ifelse(train.intent=="No","Submaximal Intended Concentric Velocity","Maximal Intended Concentric Velocity")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### rom type of strength test
secondary.str.rom.moderator.isometric.isokinetic.isotonic=secondary.str.rom.best.model.prep%>%
  update(~. + avg.rir*isometric.isokinetic.isotonic)%>%
  robust.rma.mv(cluster = secondary.data.str$study)

secondary.str.rom.moderator.isometric.isokinetic.isotonic.emms=secondary.str.rom.moderator.isometric.isokinetic.isotonic%>%emmprep(at=list(avg.rir=0:10))%>%
  emmeans(~avg.rir+isometric.isokinetic.isotonic,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.str.rom.moderator.isometric.isokinetic.isotonic)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

secondary.str.rom.moderator.isometric.isokinetic.isotonic.contrast=emmprep(secondary.str.rom.moderator.isometric.isokinetic.isotonic,
                                                                       at=list(avg.rir=secondary.str.mean.rir))%>%
  emmeans(~avg.rir+isometric.isokinetic.isotonic,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

secondary.str.rom.moderator.isometric.isokinetic.isotonic.contrast.label=paste0('Interaction Contrast (Isokinetic - Isometric) = ',
                                                                            bquote(.(round(secondary.str.rom.moderator.isometric.isokinetic.isotonic.contrast[[1,3]],5))),
                                                                            " [95% CI: ",bquote(.(round(secondary.str.rom.moderator.isometric.isokinetic.isotonic.contrast[[1,6]],5))),
                                                                            ", ",bquote(.(round(secondary.str.rom.moderator.isometric.isokinetic.isotonic.contrast[[1,7]],5))),"]\n",
                                                                            'Interaction Contrast (Isokinetic - Isotonic) = ',
                                                                            bquote(.(round(secondary.str.rom.moderator.isometric.isokinetic.isotonic.contrast[[2,3]],5))),
                                                                            " [95% CI: ",bquote(.(round(secondary.str.rom.moderator.isometric.isokinetic.isotonic.contrast[[2,6]],5))),
                                                                            ", ",bquote(.(round(secondary.str.rom.moderator.isometric.isokinetic.isotonic.contrast[[2,7]],5))),"]\n",
                                                                            'Interaction Contrast (Isometric - Isotonic) = ',
                                                                            bquote(.(round(secondary.str.rom.moderator.isometric.isokinetic.isotonic.contrast[[3,3]],5))),
                                                                            " [95% CI: ",bquote(.(round(secondary.str.rom.moderator.isometric.isokinetic.isotonic.contrast[[3,6]],5))),
                                                                            ", ",bquote(.(round(secondary.str.rom.moderator.isometric.isokinetic.isotonic.contrast[[3,7]],5))),"]")


secondary.str.rom.moderator.isometric.isokinetic.isotonic.plot=ggplot(data = secondary.data.str%>%
                                                                    mutate(isometric.isokinetic.isotonic=ifelse(isometric.isokinetic.isotonic=="isokinetic","Isokinetic",
                                                                                                                ifelse(isometric.isokinetic.isotonic=="isometric","Isometric","Isotonic"))),
                                                                  aes(x=avg.rir,
                                                                      y=yi.rom.exp))+THEME+theme(legend.position = "right",
                                                                                         plot.caption = element_text(face = "italic"))+
  facet_wrap(~isometric.isokinetic.isotonic)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Type of Strength Test x <10 Estimated RIR",
       caption = secondary.str.rom.moderator.isometric.isokinetic.isotonic.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = secondary.str.rom.moderator.isometric.isokinetic.isotonic.emms%>%
                mutate(isometric.isokinetic.isotonic=ifelse(isometric.isokinetic.isotonic=="isokinetic","Isokinetic",
                                                            ifelse(isometric.isokinetic.isotonic=="isometric","Isometric","Isotonic"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = secondary.str.rom.moderator.isometric.isokinetic.isotonic.emms%>%
                mutate(isometric.isokinetic.isotonic=ifelse(isometric.isokinetic.isotonic=="isokinetic","Isokinetic",
                                                            ifelse(isometric.isokinetic.isotonic=="isometric","Isometric","Isotonic"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = secondary.str.rom.moderator.isometric.isokinetic.isotonic.emms%>%
              mutate(isometric.isokinetic.isotonic=ifelse(isometric.isokinetic.isotonic=="isokinetic","Isokinetic",
                                                          ifelse(isometric.isokinetic.isotonic=="isometric","Isometric","Isotonic"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### rom 1RM or non-1RM
secondary.str.rom.moderator.RM.max.submax=secondary.str.rom.best.model.prep%>%
  update(~. + avg.rir*RM.max.submax)%>%
  robust.rma.mv(cluster = secondary.data.str$study)

secondary.str.rom.moderator.RM.max.submax.emms=secondary.str.rom.moderator.RM.max.submax%>%emmprep(at=list(avg.rir=0:10))%>%
  emmeans(~avg.rir+RM.max.submax,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.str.rom.moderator.RM.max.submax)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

secondary.str.rom.moderator.RM.max.submax.contrast=emmprep(secondary.str.rom.moderator.RM.max.submax,
                                                       at=list(avg.rir=secondary.str.mean.rir))%>%
  emmeans(~avg.rir+RM.max.submax,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

secondary.str.rom.moderator.RM.max.submax.contrast.label=paste0('Interaction Contrast (1RM - Non 1RM) = ',
                                                            bquote(.(round(secondary.str.rom.moderator.RM.max.submax.contrast[[1,3]],5))),
                                                            " [95% CI: ",bquote(.(round(secondary.str.rom.moderator.RM.max.submax.contrast[[1,6]],5))),
                                                            ", ",bquote(.(round(secondary.str.rom.moderator.RM.max.submax.contrast[[1,7]],5))),"]")


secondary.str.rom.moderator.RM.max.submax.plot=ggplot(data = secondary.data.str%>%
                                                    mutate(RM.max.submax=ifelse(RM.max.submax=="1RM","1RM","Non-1RM"))%>%
                                                    filter(RM.max.submax!="NA"),
                                                  aes(x=avg.rir,
                                                      y=yi.rom.exp))+THEME+theme(legend.position = "right",
                                                                         plot.caption = element_text(face = "italic"))+
  facet_wrap(~RM.max.submax)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Type of Isotonic Strength Test x <10 Estimated RIR",
       caption = secondary.str.rom.moderator.RM.max.submax.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#c1272d",
             alpha=0.75)+
  geom_ribbon(data = secondary.str.rom.moderator.RM.max.submax.emms%>%
                mutate(RM.max.submax=ifelse(RM.max.submax=="1RM","1RM","Non-1RM")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#c1272d")+
  geom_ribbon(data = secondary.str.rom.moderator.RM.max.submax.emms%>%
                mutate(RM.max.submax=ifelse(RM.max.submax=="1RM","1RM","Non-1RM")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#c1272d",alpha=0.15)+
  geom_line(data = secondary.str.rom.moderator.RM.max.submax.emms%>%
              mutate(RM.max.submax=ifelse(RM.max.submax=="1RM","1RM","Non-1RM")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#c1272d")+
  scale_size_continuous(guide="none")

#### combine rom plots
secondary.str.rom.moderator.doc=secondary.str.rom.moderator.load.plot/plot_spacer()/secondary.str.rom.moderator.volume.method.plot/plot_spacer()/
  secondary.str.rom.moderator.train.status.plot/plot_spacer()/secondary.str.rom.moderator.weeks.plot/plot_spacer()/
  secondary.str.rom.moderator.adj.sets.plot/plot_spacer()/secondary.str.rom.moderator.frequency.per.muscle.plot/plot_spacer()/
  secondary.str.rom.moderator.age.plot/plot_spacer()/secondary.str.rom.moderator.sex.percent.male.plot/plot_spacer()/
  secondary.str.rom.moderator.within.between.design.plot/plot_spacer()/secondary.str.rom.moderator.upper.lower.other.plot/plot_spacer()/
  secondary.str.rom.moderator.train.exercise.plot/plot_spacer()/secondary.str.rom.moderator.formal.cardio.intervention.plot/plot_spacer()/
  secondary.str.rom.moderator.progression.plot/plot_spacer()/secondary.str.rom.moderator.failure.definition.plot/plot_spacer()/
  secondary.str.rom.moderator.alternative.set.structure.plot/plot_spacer()/secondary.str.rom.moderator.train.intent.plot/plot_spacer()/
  secondary.str.rom.moderator.isometric.isokinetic.isotonic.plot/plot_spacer()/secondary.str.rom.moderator.RM.max.submax.plot&
  plot_annotation(title = "Moderator Analyses from Linear Fit Multilevel Model for Maximal Strength (Response Ratio)")&
  THEME

ggsave(secondary.str.rom.moderator.doc,device = "pdf",filename="secondary.str.rom.moderator.doc.pdf",width = 12,height = 156,limitsize = FALSE,
       path = "/Users/zacrobinson/Documents/For Publication PTF Meta/Plots")

### hypertrophy
#### g load
secondary.hyp.g.moderator.load=secondary.hyp.g.best.model.prep%>%
  update(~. + avg.rir:load.set)%>%
  robust.rma.mv(cluster = secondary.data.hyp$study)

secondary.hyp.g.moderator.load.emms=secondary.hyp.g.moderator.load%>%
  emmprep(at=list(avg.rir=0:10,
                  load.set=c(30,60,90)))%>%
  emmeans(~avg.rir+load.set,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.hyp.g.moderator.load)

secondary.hyp.g.moderator.load.contrast=emmprep(secondary.hyp.g.moderator.load,
                                              at=list(avg.rir=secondary.hyp.mean.rir,
                                                      load.set=secondary.hyp.mean.load.set))%>%
  emmeans(~avg.rir+load.set,
          weights = "prop")%>%
  contrast(interaction= c("consec","pairwise"))%>%
  confint()

secondary.hyp.g.moderator.load.contrast.label=paste0('Interaction Contrast (Load) = ',
                                                   bquote(.(round(secondary.hyp.g.moderator.load.contrast[[3]],5))),
                                                   " [95% CI: ",bquote(.(round(secondary.hyp.g.moderator.load.contrast[[6]],5))),
                                                   ", ",bquote(.(round(secondary.hyp.g.moderator.load.contrast[[7]],5))),"]")

secondary.hyp.g.moderator.load.plot=ggplot(data = secondary.data.hyp%>%
                                           mutate(load.set=ifelse(load.set<45,"30% of 1RM",
                                                                  ifelse(load.set>45&load.set<75,"60% of 1RM","90% of 1RM"))),
                                         aes(x=avg.rir,
                                             y=yi))+THEME+theme(plot.caption = element_text(face = "italic"))+
  facet_wrap(~load.set)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Load (% of 1RM) x <10 Estimated RIR",
       caption = secondary.hyp.g.moderator.load.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = secondary.hyp.g.moderator.load.emms%>%
                mutate(load.set=ifelse(load.set==30,"30% of 1RM",
                                       ifelse(load.set==60,"60% of 1RM","90% of 1RM"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = secondary.hyp.g.moderator.load.emms%>%
                mutate(load.set=ifelse(load.set==30,"30% of 1RM",
                                       ifelse(load.set==60,"60% of 1RM","90% of 1RM"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = secondary.hyp.g.moderator.load.emms%>%
              mutate(load.set=ifelse(load.set==30,"30% of 1RM",
                                     ifelse(load.set==60,"60% of 1RM","90% of 1RM"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### g volume equating method
secondary.hyp.g.moderator.volume.method=secondary.hyp.g.best.model.prep%>%
  update(~. + avg.rir:set.rep.equated)%>%
  robust.rma.mv(cluster = secondary.data.hyp$study)

secondary.hyp.g.moderator.volume.method.emms=secondary.hyp.g.moderator.volume.method%>%emmprep(at=list(avg.rir=0:10))%>%
  emmeans(~avg.rir+set.rep.equated,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.hyp.g.moderator.volume.method)

secondary.hyp.g.moderator.volume.method.contrast=emmprep(secondary.hyp.g.moderator.volume.method,
                                                       at=list(avg.rir=secondary.hyp.mean.rir))%>%
  emmeans(~avg.rir+set.rep.equated,
          weights = "prop")%>%
  contrast(interaction= c("consec","pairwise"))%>%
  confint()

secondary.hyp.g.moderator.volume.method.contrast.label=paste0('Interaction Contrast (Set - Rep) = ',
                                                            bquote(.(round(secondary.hyp.g.moderator.volume.method.contrast[[1,3]],5))),
                                                            " [95% CI: ",bquote(.(round(secondary.hyp.g.moderator.volume.method.contrast[[1,6]],5))),
                                                            ", ",bquote(.(round(secondary.hyp.g.moderator.volume.method.contrast[[1,7]],5))),"]\n",
                                                            'Interaction Contrast (Set - Both) = ',
                                                            bquote(.(round(secondary.hyp.g.moderator.volume.method.contrast[[2,3]],5))),
                                                            " [95% CI: ",bquote(.(round(secondary.hyp.g.moderator.volume.method.contrast[[2,6]],5))),
                                                            ", ",bquote(.(round(secondary.hyp.g.moderator.volume.method.contrast[[2,7]],5))),"]\n",
                                                            'Interaction Contrast (Rep - Both) = ',
                                                            bquote(.(round(secondary.hyp.g.moderator.volume.method.contrast[[3,3]],5))),
                                                            " [95% CI: ",bquote(.(round(secondary.hyp.g.moderator.volume.method.contrast[[3,6]],5))),
                                                            ", ",bquote(.(round(secondary.hyp.g.moderator.volume.method.contrast[[3,7]],5))),"]")


secondary.hyp.g.moderator.volume.method.plot=ggplot(data = secondary.data.hyp%>%
                                                    mutate(set.rep.equated=ifelse(set.rep.equated=="set","Set Equated",
                                                                                  ifelse(set.rep.equated=="rep","Repetition Equated","Both Equated"))),
                                                  aes(x=avg.rir,
                                                      y=yi))+THEME+theme(legend.position = "right",
                                                                         plot.caption = element_text(face = "italic"))+
  facet_wrap(~set.rep.equated)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Method of Volume Equating x <10 Estimated RIR",
       caption = secondary.hyp.g.moderator.volume.method.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = secondary.hyp.g.moderator.volume.method.emms%>%
                mutate(set.rep.equated=ifelse(set.rep.equated=="set","Set Equated",
                                              ifelse(set.rep.equated=="rep","Repetition Equated","Both Equated"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = secondary.hyp.g.moderator.volume.method.emms%>%
                mutate(set.rep.equated=ifelse(set.rep.equated=="set","Set Equated",
                                              ifelse(set.rep.equated=="rep","Repetition Equated","Both Equated"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = secondary.hyp.g.moderator.volume.method.emms%>%
              mutate(set.rep.equated=ifelse(set.rep.equated=="set","Set Equated",
                                            ifelse(set.rep.equated=="rep","Repetition Equated","Both Equated"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### g training status
secondary.hyp.g.moderator.train.status=secondary.hyp.g.best.model.prep%>%
  update(~. + avg.rir:train.status)%>%
  robust.rma.mv(cluster = secondary.data.hyp$study)

secondary.hyp.g.moderator.train.status.emms=secondary.hyp.g.moderator.train.status%>%emmprep(at=list(avg.rir=0:10))%>%
  emmeans(~avg.rir+train.status,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.hyp.g.moderator.train.status)

secondary.hyp.g.moderator.train.status.contrast=emmprep(secondary.hyp.g.moderator.train.status,
                                                      at=list(avg.rir=secondary.hyp.mean.rir))%>%
  emmeans(~avg.rir+train.status,
          weights = "prop")%>%
  contrast(interaction= c("consec","pairwise"))%>%
  confint()

secondary.hyp.g.moderator.train.status.contrast.label=paste0('Interaction Contrast (Trained - Untrained) = ',
                                                           bquote(.(round(secondary.hyp.g.moderator.train.status.contrast[[1,3]],5))),
                                                           " [95% CI: ",bquote(.(round(secondary.hyp.g.moderator.train.status.contrast[[1,6]],5))),
                                                           ", ",bquote(.(round(secondary.hyp.g.moderator.train.status.contrast[[1,7]],5))),"]")


secondary.hyp.g.moderator.train.status.plot=ggplot(data = secondary.data.hyp%>%
                                                   mutate(train.status=ifelse(train.status=="trained","Trained","Untrained")),
                                                 aes(x=avg.rir,
                                                     y=yi))+THEME+theme(legend.position = "right",
                                                                        plot.caption = element_text(face = "italic"))+
  facet_wrap(~train.status)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Training Status x <10 Estimated RIR",
       caption = secondary.hyp.g.moderator.train.status.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = secondary.hyp.g.moderator.train.status.emms%>%
                mutate(train.status=ifelse(train.status=="trained","Trained","Untrained")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = secondary.hyp.g.moderator.train.status.emms%>%
                mutate(train.status=ifelse(train.status=="trained","Trained","Untrained")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = secondary.hyp.g.moderator.train.status.emms%>%
              mutate(train.status=ifelse(train.status=="trained","Trained","Untrained")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### g weeks
secondary.hyp.g.moderator.weeks=secondary.hyp.g.best.model.prep%>%
  update(~. + avg.rir:weeks)%>%
  robust.rma.mv(cluster = secondary.data.hyp$study)

secondary.hyp.g.moderator.weeks.emms=secondary.hyp.g.moderator.weeks%>%
  emmprep(at=list(avg.rir=0:10,
                  weeks=c(4,8,12)))%>%
  emmeans(~avg.rir+weeks,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.hyp.g.moderator.weeks)

secondary.hyp.g.moderator.weeks.contrast=emmprep(secondary.hyp.g.moderator.weeks,
                                               at=list(avg.rir=secondary.hyp.mean.rir,
                                                       weeks=secondary.hyp.mean.weeks))%>%
  emmeans(~avg.rir+weeks,
          weights = "prop")%>%
  contrast("consec",interaction=TRUE)%>%
  confint()

secondary.hyp.g.moderator.weeks.contrast.label=paste0('Interaction Contrast (Weeks) = ',
                                                    bquote(.(round(secondary.hyp.g.moderator.weeks.contrast[[3]],5))),
                                                    " [95% CI: ",bquote(.(round(secondary.hyp.g.moderator.weeks.contrast[[6]],5))),
                                                    ", ",bquote(.(round(secondary.hyp.g.moderator.weeks.contrast[[7]],5))),"]")

secondary.hyp.g.moderator.weeks.plot=ggplot(data = secondary.data.hyp%>%
                                            mutate(weeks=ifelse(weeks<6,"4 Weeks",
                                                                ifelse(weeks>6&weeks<10,"8 Weeks","12 Weeks"))),
                                          aes(x=avg.rir,
                                              y=yi))+THEME+theme(plot.caption = element_text(face = "italic"))+
  facet_wrap(~fct_relevel(weeks,c("4 Weeks","8 Weeks","12 Weeks")))+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Intervention Duration (Weeks) x <10 Estimated RIR",
       caption = secondary.hyp.g.moderator.weeks.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = secondary.hyp.g.moderator.weeks.emms%>%
                mutate(weeks=ifelse(weeks==4,"4 Weeks",
                                    ifelse(weeks==8,"8 Weeks","12 Weeks"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = secondary.hyp.g.moderator.weeks.emms%>%
                mutate(weeks=ifelse(weeks==4,"4 Weeks",
                                    ifelse(weeks==8,"8 Weeks","12 Weeks"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = secondary.hyp.g.moderator.weeks.emms%>%
              mutate(weeks=ifelse(weeks==4,"4 Weeks",
                                  ifelse(weeks==8,"8 Weeks","12 Weeks"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### g adjusted sets 
secondary.hyp.g.moderator.adj.sets=secondary.hyp.g.best.model.prep%>%
  update(~. + avg.rir*adj.sets.week)%>%
  robust.rma.mv(cluster = secondary.data.hyp$study)

secondary.hyp.g.moderator.adj.sets.emms=secondary.hyp.g.moderator.adj.sets%>%
  emmprep(at=list(avg.rir=0:10,
                  adj.sets.week=c(8,16,24)))%>%
  emmeans(~avg.rir+adj.sets.week,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.hyp.g.moderator.adj.sets)

secondary.hyp.g.moderator.adj.sets.contrast=emmprep(secondary.hyp.g.moderator.adj.sets,
                                                  at=list(avg.rir=secondary.hyp.mean.rir,
                                                          adj.sets.week=secondary.hyp.mean.adj.sets.week))%>%
  emmeans(~avg.rir+adj.sets.week,
          weights = "prop")%>%
  contrast("consec",interaction=TRUE)%>%
  confint()

secondary.hyp.g.moderator.adj.sets.contrast.label=paste0('Interaction Contrast (Adjusted Direct Sets) = ',
                                                       bquote(.(round(secondary.hyp.g.moderator.adj.sets.contrast[[3]],5))),
                                                       " [95% CI: ",bquote(.(round(secondary.hyp.g.moderator.adj.sets.contrast[[6]],5))),
                                                       ", ",bquote(.(round(secondary.hyp.g.moderator.adj.sets.contrast[[7]],5))),"]")

secondary.hyp.g.moderator.adj.sets.plot=ggplot(data = secondary.data.hyp%>%
                                               mutate(adj.sets.week=ifelse(adj.sets.week<12,"8 Sets",
                                                                           ifelse(adj.sets.week>12&adj.sets.week<20,"16 Sets","24 Sets"))),
                                             aes(x=avg.rir,
                                                 y=yi))+THEME+theme(plot.caption = element_text(face = "italic"))+
  facet_wrap(~fct_relevel(adj.sets.week,c("8 Sets","16 Sets","24 Sets")))+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Adjusted Direct Sets Per Week Per Muscle Group x <10 Estimated RIR",
       caption = secondary.hyp.g.moderator.adj.sets.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = secondary.hyp.g.moderator.adj.sets.emms%>%
                mutate(adj.sets.week=ifelse(adj.sets.week==8,"8 Sets",
                                            ifelse(adj.sets.week==16,"16 Sets","24 Sets"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = secondary.hyp.g.moderator.adj.sets.emms%>%
                mutate(adj.sets.week=ifelse(adj.sets.week==8,"8 Sets",
                                            ifelse(adj.sets.week==16,"16 Sets","24 Sets"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = secondary.hyp.g.moderator.adj.sets.emms%>%
              mutate(adj.sets.week=ifelse(adj.sets.week==8,"8 Sets",
                                          ifelse(adj.sets.week==16,"16 Sets","24 Sets"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### g frequency
secondary.hyp.g.moderator.frequency.per.muscle=secondary.hyp.g.best.model.prep%>%
  update(~. + avg.rir*frequency.per.muscle)%>%
  robust.rma.mv(cluster = secondary.data.hyp$study)

secondary.hyp.g.moderator.frequency.per.muscle.emms=secondary.hyp.g.moderator.frequency.per.muscle%>%
  emmprep(at=list(avg.rir=0:10,
                  frequency.per.muscle=c(1,2,3)))%>%
  emmeans(~avg.rir+frequency.per.muscle,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.hyp.g.moderator.frequency.per.muscle)

secondary.hyp.g.moderator.frequency.per.muscle.contrast=emmprep(secondary.hyp.g.moderator.frequency.per.muscle,
                                                              at=list(avg.rir=secondary.hyp.mean.rir,
                                                                      frequency.per.muscle=secondary.hyp.mean.frequency.per.muscle))%>%
  emmeans(~avg.rir+frequency.per.muscle,
          weights = "prop")%>%
  contrast("consec",interaction=TRUE)%>%
  confint()

secondary.hyp.g.moderator.frequency.per.muscle.contrast.label=paste0('Interaction Contrast (Direct Frequency) = ',
                                                                   bquote(.(round(secondary.hyp.g.moderator.frequency.per.muscle.contrast[[3]],5))),
                                                                   " [95% CI: ",bquote(.(round(secondary.hyp.g.moderator.frequency.per.muscle.contrast[[6]],5))),
                                                                   ", ",bquote(.(round(secondary.hyp.g.moderator.frequency.per.muscle.contrast[[7]],5))),"]")

secondary.hyp.g.moderator.frequency.per.muscle.plot=ggplot(data = secondary.data.hyp%>%
                                                           mutate(frequency.per.muscle=ifelse(frequency.per.muscle<1.5,"1x Per Week",
                                                                                              ifelse(frequency.per.muscle>1.5&frequency.per.muscle<2.5,"2x Per Week","3x Per Week"))),
                                                         aes(x=avg.rir,
                                                             y=yi))+THEME+theme(plot.caption = element_text(face = "italic"))+
  facet_wrap(~frequency.per.muscle)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Direct Frequency Per Week Per Muscle Group x <10 Estimated RIR",
       caption = secondary.hyp.g.moderator.frequency.per.muscle.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = secondary.hyp.g.moderator.frequency.per.muscle.emms%>%
                mutate(frequency.per.muscle=ifelse(frequency.per.muscle==1,"1x Per Week",
                                                   ifelse(frequency.per.muscle==2,"2x Per Week","3x Per Week"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = secondary.hyp.g.moderator.frequency.per.muscle.emms%>%
                mutate(frequency.per.muscle=ifelse(frequency.per.muscle==1,"1x Per Week",
                                                   ifelse(frequency.per.muscle==2,"2x Per Week","3x Per Week"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = secondary.hyp.g.moderator.frequency.per.muscle.emms%>%
              mutate(frequency.per.muscle=ifelse(frequency.per.muscle==1,"1x Per Week",
                                                 ifelse(frequency.per.muscle==2,"2x Per Week","3x Per Week"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### g age
secondary.hyp.g.moderator.age=secondary.hyp.g.best.model.prep%>%
  update(~. + avg.rir*age)%>%
  robust.rma.mv(cluster = secondary.data.hyp$study)

secondary.hyp.g.moderator.age.emms=secondary.hyp.g.moderator.age%>%
  emmprep(at=list(avg.rir=0:10,
                  age=c(20,40,60)))%>%
  emmeans(~avg.rir+age,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.hyp.g.moderator.age)

secondary.hyp.g.moderator.age.contrast=emmprep(secondary.hyp.g.moderator.age,
                                             at=list(avg.rir=secondary.hyp.mean.rir,
                                                     age=secondary.hyp.mean.age))%>%
  emmeans(~avg.rir+age,
          weights = "prop")%>%
  contrast("consec",interaction=TRUE)%>%
  confint()

secondary.hyp.g.moderator.age.contrast.label=paste0('Interaction Contrast (Age) = ',
                                                  bquote(.(round(secondary.hyp.g.moderator.age.contrast[[3]],5))),
                                                  " [95% CI: ",bquote(.(round(secondary.hyp.g.moderator.age.contrast[[6]],5))),
                                                  ", ",bquote(.(round(secondary.hyp.g.moderator.age.contrast[[7]],5))),"]")

secondary.hyp.g.moderator.age.plot=ggplot(data = secondary.data.hyp%>%
                                          mutate(age=ifelse(age<30,"20 Years Old",
                                                            ifelse(age>30&age<50,"40 Years Old","60 Years Old"))),
                                        aes(x=avg.rir,
                                            y=yi))+THEME+theme(plot.caption = element_text(face = "italic"))+
  facet_wrap(~age)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Age x <10 Estimated RIR",
       caption = secondary.hyp.g.moderator.age.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = secondary.hyp.g.moderator.age.emms%>%
                mutate(age=ifelse(age==20,"20 Years Old",
                                  ifelse(age==40,"40 Years Old","60 Years Old"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = secondary.hyp.g.moderator.age.emms%>%
                mutate(age=ifelse(age==20,"20 Years Old",
                                  ifelse(age==40,"40 Years Old","60 Years Old"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = secondary.hyp.g.moderator.age.emms%>%
              mutate(age=ifelse(age==20,"20 Years Old",
                                ifelse(age==40,"40 Years Old","60 Years Old"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### g sex
secondary.hyp.g.moderator.sex.percent.male=secondary.hyp.g.best.model.prep%>%
  update(~. + avg.rir*sex.percent.male)%>%
  robust.rma.mv(cluster = secondary.data.hyp$study)

secondary.hyp.g.moderator.sex.percent.male.emms=secondary.hyp.g.moderator.sex.percent.male%>%
  emmprep(at=list(avg.rir=0:10,
                  sex.percent.male=c(0,100)))%>%
  emmeans(~avg.rir+sex.percent.male,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.hyp.g.moderator.sex.percent.male)

secondary.hyp.g.moderator.sex.percent.male.contrast=emmprep(secondary.hyp.g.moderator.sex.percent.male,
                                                          at=list(avg.rir=secondary.hyp.mean.rir,
                                                                  sex.percent.male=secondary.hyp.mean.sex))%>%
  emmeans(~avg.rir+sex.percent.male,
          weights = "prop")%>%
  contrast("consec",interaction=TRUE)%>%
  confint()

secondary.hyp.g.moderator.sex.percent.male.contrast.label=paste0('Interaction Contrast (Sex) = ',
                                                               bquote(.(round(secondary.hyp.g.moderator.sex.percent.male.contrast[[3]],5))),
                                                               " [95% CI: ",bquote(.(round(secondary.hyp.g.moderator.sex.percent.male.contrast[[6]],5))),
                                                               ", ",bquote(.(round(secondary.hyp.g.moderator.sex.percent.male.contrast[[7]],5))),"]")

secondary.hyp.g.moderator.sex.percent.male.plot=ggplot(data = secondary.data.hyp%>%
                                                       mutate(sex.percent.male=ifelse(sex.percent.male<50,"Female","Male")),
                                                     aes(x=avg.rir,
                                                         y=yi))+THEME+theme(plot.caption = element_text(face = "italic"))+
  facet_wrap(~sex.percent.male)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Sex (% Male) x <10 Estimated RIR",
       caption = secondary.hyp.g.moderator.sex.percent.male.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = secondary.hyp.g.moderator.sex.percent.male.emms%>%
                mutate(sex.percent.male=ifelse(sex.percent.male==0,"Female","Male")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = secondary.hyp.g.moderator.sex.percent.male.emms%>%
                mutate(sex.percent.male=ifelse(sex.percent.male==0,"Female","Male")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = secondary.hyp.g.moderator.sex.percent.male.emms%>%
              mutate(sex.percent.male=ifelse(sex.percent.male==0,"Female","Male")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### g study design
secondary.hyp.g.moderator.within.between.design=secondary.hyp.g.best.model.prep%>%
  update(~. + avg.rir*within.between.design)%>%
  robust.rma.mv(cluster = secondary.data.hyp$study)

secondary.hyp.g.moderator.within.between.design.emms=secondary.hyp.g.moderator.within.between.design%>%emmprep(at=list(avg.rir=0:10))%>%
  emmeans(~avg.rir+within.between.design,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.hyp.g.moderator.within.between.design)

secondary.hyp.g.moderator.within.between.design.contrast=emmprep(secondary.hyp.g.moderator.within.between.design,
                                                               at=list(avg.rir=secondary.hyp.mean.rir))%>%
  emmeans(~avg.rir+within.between.design,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()

secondary.hyp.g.moderator.within.between.design.contrast.label=paste0('Interaction Contrast (Between - Within) = ',
                                                                    bquote(.(round(secondary.hyp.g.moderator.within.between.design.contrast[[1,3]],5))),
                                                                    " [95% CI: ",bquote(.(round(secondary.hyp.g.moderator.within.between.design.contrast[[1,6]],5))),
                                                                    ", ",bquote(.(round(secondary.hyp.g.moderator.within.between.design.contrast[[1,7]],5))),"]")


secondary.hyp.g.moderator.within.between.design.plot=ggplot(data = secondary.data.hyp%>%
                                                            mutate(within.between.design=ifelse(within.between.design=="Between","Between Subjects","Within Subjects")),
                                                          aes(x=avg.rir,
                                                              y=yi))+THEME+theme(legend.position = "right",
                                                                                 plot.caption = element_text(face = "italic"))+
  facet_wrap(~within.between.design)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Study Design x <10 Estimated RIR",
       caption = secondary.hyp.g.moderator.within.between.design.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = secondary.hyp.g.moderator.within.between.design.emms%>%
                mutate(within.between.design=ifelse(within.between.design=="Between","Between Subjects","Within Subjects")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = secondary.hyp.g.moderator.within.between.design.emms%>%
                mutate(within.between.design=ifelse(within.between.design=="Between","Between Subjects","Within Subjects")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = secondary.hyp.g.moderator.within.between.design.emms%>%
              mutate(within.between.design=ifelse(within.between.design=="Between","Between Subjects","Within Subjects")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### g upper / lower outcome
secondary.hyp.g.moderator.upper.lower.other=secondary.hyp.g.best.model.prep%>%
  update(~. + avg.rir*upper.lower.other)%>%
  robust.rma.mv(cluster = secondary.data.hyp$study)

secondary.hyp.g.moderator.upper.lower.other.emms=secondary.hyp.g.moderator.upper.lower.other%>%emmprep(at=list(avg.rir=0:10))%>%
  emmeans(~avg.rir+upper.lower.other,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.hyp.g.moderator.upper.lower.other)

secondary.hyp.g.moderator.upper.lower.other.contrast=emmprep(secondary.hyp.g.moderator.upper.lower.other,
                                                           at=list(avg.rir=secondary.hyp.mean.rir))%>%
  emmeans(~avg.rir+upper.lower.other,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()

secondary.hyp.g.moderator.upper.lower.other.contrast.label=paste0('Interaction Contrast (Lower - Upper) = ',
                                                                bquote(.(round(secondary.hyp.g.moderator.upper.lower.other.contrast[[1,3]],5))),
                                                                " [95% CI: ",bquote(.(round(secondary.hyp.g.moderator.upper.lower.other.contrast[[1,6]],5))),
                                                                ", ",bquote(.(round(secondary.hyp.g.moderator.upper.lower.other.contrast[[1,7]],5))),"]")


secondary.hyp.g.moderator.upper.lower.other.plot=ggplot(data = secondary.data.hyp%>%
                                                        mutate(upper.lower.other=ifelse(upper.lower.other=="Lower","Lower Body Outcomes","Upper Body Outcomes")),
                                                      aes(x=avg.rir,
                                                          y=yi))+THEME+theme(legend.position = "right",
                                                                             plot.caption = element_text(face = "italic"))+
  facet_wrap(~upper.lower.other)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Upper/Lower Body Outcome x <10 Estimated RIR",
       caption = secondary.hyp.g.moderator.upper.lower.other.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = secondary.hyp.g.moderator.upper.lower.other.emms%>%
                mutate(upper.lower.other=ifelse(upper.lower.other=="Lower","Lower Body Outcomes","Upper Body Outcomes")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = secondary.hyp.g.moderator.upper.lower.other.emms%>%
                mutate(upper.lower.other=ifelse(upper.lower.other=="Lower","Lower Body Outcomes","Upper Body Outcomes")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = secondary.hyp.g.moderator.upper.lower.other.emms%>%
              mutate(upper.lower.other=ifelse(upper.lower.other=="Lower","Lower Body Outcomes","Upper Body Outcomes")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### g exercise selection
secondary.hyp.g.moderator.train.exercise=secondary.hyp.g.best.model.prep%>%
  update(~. + avg.rir*train.exercise)%>%
  robust.rma.mv(cluster = secondary.data.hyp$study)

secondary.hyp.g.moderator.train.exercise.emms=secondary.hyp.g.moderator.train.exercise%>%emmprep(at=list(avg.rir=0:10))%>%
  emmeans(~avg.rir+train.exercise,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.hyp.g.moderator.train.exercise)

secondary.hyp.g.moderator.train.exercise.contrast=emmprep(secondary.hyp.g.moderator.train.exercise,
                                                        at=list(avg.rir=secondary.hyp.mean.rir))%>%
  emmeans(~avg.rir+train.exercise,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()

secondary.hyp.g.moderator.train.exercise.contrast.label=paste0('Interaction Contrast (Both - Multi-joint) = ',
                                                             bquote(.(round(secondary.hyp.g.moderator.train.exercise.contrast[[1,3]],5))),
                                                             " [95% CI: ",bquote(.(round(secondary.hyp.g.moderator.train.exercise.contrast[[1,6]],5))),
                                                             ", ",bquote(.(round(secondary.hyp.g.moderator.train.exercise.contrast[[1,7]],5))),"]\n",
                                                             'Interaction Contrast (Both - Single-joint) = ',
                                                             bquote(.(round(secondary.hyp.g.moderator.train.exercise.contrast[[2,3]],5))),
                                                             " [95% CI: ",bquote(.(round(secondary.hyp.g.moderator.train.exercise.contrast[[2,6]],5))),
                                                             ", ",bquote(.(round(secondary.hyp.g.moderator.train.exercise.contrast[[2,7]],5))),"]\n",
                                                             'Interaction Contrast (Multi-joint - Single-joint) = ',
                                                             bquote(.(round(secondary.hyp.g.moderator.train.exercise.contrast[[3,3]],5))),
                                                             " [95% CI: ",bquote(.(round(secondary.hyp.g.moderator.train.exercise.contrast[[3,6]],5))),
                                                             ", ",bquote(.(round(secondary.hyp.g.moderator.train.exercise.contrast[[3,7]],5))),"]")


secondary.hyp.g.moderator.train.exercise.plot=ggplot(data = secondary.data.hyp%>%
                                                     mutate(train.exercise=ifelse(train.exercise=="Both","Both",
                                                                                  ifelse(train.exercise=="Multi-joint","Multi-joint Exercises","Single-joint Exercises"))),
                                                   aes(x=avg.rir,
                                                       y=yi))+THEME+theme(legend.position = "right",
                                                                          plot.caption = element_text(face = "italic"))+
  facet_wrap(~train.exercise)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Exercise Selection x <10 Estimated RIR",
       caption = secondary.hyp.g.moderator.train.exercise.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = secondary.hyp.g.moderator.train.exercise.emms%>%
                mutate(train.exercise=ifelse(train.exercise=="Both","Both",
                                             ifelse(train.exercise=="Multi-joint","Multi-joint Exercises","Single-joint Exercises"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = secondary.hyp.g.moderator.train.exercise.emms%>%
                mutate(train.exercise=ifelse(train.exercise=="Both","Both",
                                             ifelse(train.exercise=="Multi-joint","Multi-joint Exercises","Single-joint Exercises"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = secondary.hyp.g.moderator.train.exercise.emms%>%
              mutate(train.exercise=ifelse(train.exercise=="Both","Both",
                                           ifelse(train.exercise=="Multi-joint","Multi-joint Exercises","Single-joint Exercises"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### g concurrent training
secondary.hyp.g.moderator.formal.cardio.intervention=secondary.hyp.g.best.model.prep%>%
  update(~. + avg.rir*formal.cardio.intervention)%>%
  robust.rma.mv(cluster = secondary.data.hyp$study)

secondary.hyp.g.moderator.formal.cardio.intervention.emms=secondary.hyp.g.moderator.formal.cardio.intervention%>%emmprep(at=list(avg.rir=0:10))%>%
  emmeans(~avg.rir+formal.cardio.intervention,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.hyp.g.moderator.formal.cardio.intervention)

secondary.hyp.g.moderator.formal.cardio.intervention.contrast=emmprep(secondary.hyp.g.moderator.formal.cardio.intervention,
                                                                    at=list(avg.rir=secondary.hyp.mean.rir))%>%
  emmeans(~avg.rir+formal.cardio.intervention,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()

secondary.hyp.g.moderator.formal.cardio.intervention.contrast.label=paste0('Interaction Contrast (RT Only - Concurrent) = ',
                                                                         bquote(.(round(secondary.hyp.g.moderator.formal.cardio.intervention.contrast[[1,3]],5))),
                                                                         " [95% CI: ",bquote(.(round(secondary.hyp.g.moderator.formal.cardio.intervention.contrast[[1,6]],5))),
                                                                         ", ",bquote(.(round(secondary.hyp.g.moderator.formal.cardio.intervention.contrast[[1,7]],5))),"]")


secondary.hyp.g.moderator.formal.cardio.intervention.plot=ggplot(data = secondary.data.hyp%>%
                                                                 mutate(formal.cardio.intervention=ifelse(formal.cardio.intervention=="No","Only Resistance Training","Concurrent Training")),
                                                               aes(x=avg.rir,
                                                                   y=yi))+THEME+theme(legend.position = "right",
                                                                                      plot.caption = element_text(face = "italic"))+
  facet_wrap(~formal.cardio.intervention)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Concurrent Training x <10 Estimated RIR",
       caption = secondary.hyp.g.moderator.formal.cardio.intervention.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = secondary.hyp.g.moderator.formal.cardio.intervention.emms%>%
                mutate(formal.cardio.intervention=ifelse(formal.cardio.intervention=="No","Only Resistance Training","Concurrent Training")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = secondary.hyp.g.moderator.formal.cardio.intervention.emms%>%
                mutate(formal.cardio.intervention=ifelse(formal.cardio.intervention=="No","Only Resistance Training","Concurrent Training")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = secondary.hyp.g.moderator.formal.cardio.intervention.emms%>%
              mutate(formal.cardio.intervention=ifelse(formal.cardio.intervention=="No","Only Resistance Training","Concurrent Training")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### g progressive overload
secondary.hyp.g.moderator.progression=secondary.hyp.g.best.model.prep%>%
  update(~. + avg.rir*progression)%>%
  robust.rma.mv(cluster = secondary.data.hyp$study)

secondary.hyp.g.moderator.progression.emms=secondary.hyp.g.moderator.progression%>%emmprep(at=list(avg.rir=0:10))%>%
  emmeans(~avg.rir+progression,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.hyp.g.moderator.progression)

secondary.hyp.g.moderator.progression.contrast=emmprep(secondary.hyp.g.moderator.progression,
                                                     at=list(avg.rir=secondary.hyp.mean.rir))%>%
  emmeans(~avg.rir+progression,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()

secondary.hyp.g.moderator.progression.contrast.label=paste0('Interaction Contrast (No Overload - Progressive Overload) = ',
                                                          bquote(.(round(secondary.hyp.g.moderator.progression.contrast[[1,3]],5))),
                                                          " [95% CI: ",bquote(.(round(secondary.hyp.g.moderator.progression.contrast[[1,6]],5))),
                                                          ", ",bquote(.(round(secondary.hyp.g.moderator.progression.contrast[[1,7]],5))),"]")


secondary.hyp.g.moderator.progression.plot=ggplot(data = secondary.data.hyp%>%
                                                  mutate(progression=ifelse(progression=="No","No Progressive Overload","Progressive Overload")),
                                                aes(x=avg.rir,
                                                    y=yi))+THEME+theme(legend.position = "right",
                                                                       plot.caption = element_text(face = "italic"))+
  facet_wrap(~progression)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Progressive Overload x <10 Estimated RIR",
       caption = secondary.hyp.g.moderator.progression.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = secondary.hyp.g.moderator.progression.emms%>%
                mutate(progression=ifelse(progression=="No","No Progressive Overload","Progressive Overload")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = secondary.hyp.g.moderator.progression.emms%>%
                mutate(progression=ifelse(progression=="No","No Progressive Overload","Progressive Overload")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = secondary.hyp.g.moderator.progression.emms%>%
              mutate(progression=ifelse(progression=="No","No Progressive Overload","Progressive Overload")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### g failure definition
secondary.hyp.g.moderator.failure.definition=secondary.hyp.g.best.model.prep%>%
  update(~. + avg.rir*failure.definition)%>%
  robust.rma.mv(cluster = secondary.data.hyp$study)

secondary.hyp.g.moderator.failure.definition.emms=secondary.hyp.g.moderator.failure.definition%>%emmprep(at=list(avg.rir=0:10))%>%
  emmeans(~avg.rir+failure.definition,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.hyp.g.moderator.failure.definition)

secondary.hyp.g.moderator.failure.definition.contrast=emmprep(secondary.hyp.g.moderator.failure.definition,
                                                            at=list(avg.rir=secondary.hyp.mean.rir))%>%
  emmeans(~avg.rir+failure.definition,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()

secondary.hyp.g.moderator.failure.definition.contrast.label=paste0('Interaction Contrast (No Definition - Definition Provided) = ',
                                                                 bquote(.(round(secondary.hyp.g.moderator.failure.definition.contrast[[1,3]],5))),
                                                                 " [95% CI: ",bquote(.(round(secondary.hyp.g.moderator.failure.definition.contrast[[1,6]],5))),
                                                                 ", ",bquote(.(round(secondary.hyp.g.moderator.failure.definition.contrast[[1,7]],5))),"]")


secondary.hyp.g.moderator.failure.definition.plot=ggplot(data = secondary.data.hyp%>%
                                                         mutate(failure.definition=ifelse(failure.definition=="No","No Failure Definition","Failure Definition Provided")),
                                                       aes(x=avg.rir,
                                                           y=yi))+THEME+theme(legend.position = "right",
                                                                              plot.caption = element_text(face = "italic"))+
  facet_wrap(~failure.definition)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Failure Definition x <10 Estimated RIR",
       caption = secondary.hyp.g.moderator.failure.definition.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = secondary.hyp.g.moderator.failure.definition.emms%>%
                mutate(failure.definition=ifelse(failure.definition=="No","No Failure Definition","Failure Definition Provided")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = secondary.hyp.g.moderator.failure.definition.emms%>%
                mutate(failure.definition=ifelse(failure.definition=="No","No Failure Definition","Failure Definition Provided")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = secondary.hyp.g.moderator.failure.definition.emms%>%
              mutate(failure.definition=ifelse(failure.definition=="No","No Failure Definition","Failure Definition Provided")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### g set structure
secondary.hyp.g.moderator.alternative.set.structure=secondary.hyp.g.best.model.prep%>%
  update(~. + avg.rir*alternative.set.structure)%>%
  robust.rma.mv(cluster = secondary.data.hyp$study)

secondary.hyp.g.moderator.alternative.set.structure.emms=secondary.hyp.g.moderator.alternative.set.structure%>%emmprep(at=list(avg.rir=0:10))%>%
  emmeans(~avg.rir+alternative.set.structure,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.hyp.g.moderator.alternative.set.structure)

secondary.hyp.g.moderator.alternative.set.structure.contrast=emmprep(secondary.hyp.g.moderator.alternative.set.structure,
                                                                   at=list(avg.rir=secondary.hyp.mean.rir))%>%
  emmeans(~avg.rir+alternative.set.structure,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()

secondary.hyp.g.moderator.alternative.set.structure.contrast.label=paste0('Interaction Contrast (Traditional - Alternative) = ',
                                                                        bquote(.(round(secondary.hyp.g.moderator.alternative.set.structure.contrast[[1,3]],5))),
                                                                        " [95% CI: ",bquote(.(round(secondary.hyp.g.moderator.alternative.set.structure.contrast[[1,6]],5))),
                                                                        ", ",bquote(.(round(secondary.hyp.g.moderator.alternative.set.structure.contrast[[1,7]],5))),"]")


secondary.hyp.g.moderator.alternative.set.structure.plot=ggplot(data = secondary.data.hyp%>%
                                                                mutate(alternative.set.structure=ifelse(alternative.set.structure=="No","Traditional Set Structure","Alternative Set Structure")),
                                                              aes(x=avg.rir,
                                                                  y=yi))+THEME+theme(legend.position = "right",
                                                                                     plot.caption = element_text(face = "italic"))+
  facet_wrap(~alternative.set.structure)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Set Structure x <10 Estimated RIR",
       caption = secondary.hyp.g.moderator.alternative.set.structure.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = secondary.hyp.g.moderator.alternative.set.structure.emms%>%
                mutate(alternative.set.structure=ifelse(alternative.set.structure=="No","Traditional Set Structure","Alternative Set Structure")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = secondary.hyp.g.moderator.alternative.set.structure.emms%>%
                mutate(alternative.set.structure=ifelse(alternative.set.structure=="No","Traditional Set Structure","Alternative Set Structure")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = secondary.hyp.g.moderator.alternative.set.structure.emms%>%
              mutate(alternative.set.structure=ifelse(alternative.set.structure=="No","Traditional Set Structure","Alternative Set Structure")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### g intended velocity
secondary.hyp.g.moderator.train.intent=secondary.hyp.g.best.model.prep%>%
  update(~. + avg.rir*train.intent)%>%
  robust.rma.mv(cluster = secondary.data.hyp$study)

secondary.hyp.g.moderator.train.intent.emms=secondary.hyp.g.moderator.train.intent%>%emmprep(at=list(avg.rir=0:10))%>%
  emmeans(~avg.rir+train.intent,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.hyp.g.moderator.train.intent)

secondary.hyp.g.moderator.train.intent.contrast=emmprep(secondary.hyp.g.moderator.train.intent,
                                                      at=list(avg.rir=secondary.hyp.mean.rir))%>%
  emmeans(~avg.rir+train.intent,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()

secondary.hyp.g.moderator.train.intent.contrast.label=paste0('Interaction Contrast (Submaximal - Maximal) = ',
                                                           bquote(.(round(secondary.hyp.g.moderator.train.intent.contrast[[1,3]],5))),
                                                           " [95% CI: ",bquote(.(round(secondary.hyp.g.moderator.train.intent.contrast[[1,6]],5))),
                                                           ", ",bquote(.(round(secondary.hyp.g.moderator.train.intent.contrast[[1,7]],5))),"]")


secondary.hyp.g.moderator.train.intent.plot=ggplot(data = secondary.data.hyp%>%
                                                   mutate(train.intent=ifelse(train.intent=="No","Submaximal Intended Concentric Velocity","Maximal Intended Concentric Velocity")),
                                                 aes(x=avg.rir,
                                                     y=yi))+THEME+theme(legend.position = "right",
                                                                        plot.caption = element_text(face = "italic"))+
  facet_wrap(~train.intent)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Standardized Mean Change in Strength (Hedge's g)",
       title = "Intended Concentric Velocity x <10 Estimated RIR",
       caption = secondary.hyp.g.moderator.train.intent.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = secondary.hyp.g.moderator.train.intent.emms%>%
                mutate(train.intent=ifelse(train.intent=="No","Submaximal Intended Concentric Velocity","Maximal Intended Concentric Velocity")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = secondary.hyp.g.moderator.train.intent.emms%>%
                mutate(train.intent=ifelse(train.intent=="No","Submaximal Intended Concentric Velocity","Maximal Intended Concentric Velocity")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = secondary.hyp.g.moderator.train.intent.emms%>%
              mutate(train.intent=ifelse(train.intent=="No","Submaximal Intended Concentric Velocity","Maximal Intended Concentric Velocity")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### combine g plots
secondary.hyp.g.moderator.doc=secondary.hyp.g.moderator.load.plot/plot_spacer()/secondary.hyp.g.moderator.volume.method.plot/plot_spacer()/
  secondary.hyp.g.moderator.train.status.plot/plot_spacer()/secondary.hyp.g.moderator.weeks.plot/plot_spacer()/
  secondary.hyp.g.moderator.adj.sets.plot/plot_spacer()/secondary.hyp.g.moderator.frequency.per.muscle.plot/plot_spacer()/
  secondary.hyp.g.moderator.age.plot/plot_spacer()/secondary.hyp.g.moderator.sex.percent.male.plot/plot_spacer()/
  secondary.hyp.g.moderator.within.between.design.plot/plot_spacer()/secondary.hyp.g.moderator.upper.lower.other.plot/plot_spacer()/
  secondary.hyp.g.moderator.train.exercise.plot/plot_spacer()/secondary.hyp.g.moderator.formal.cardio.intervention.plot/plot_spacer()/
  secondary.hyp.g.moderator.progression.plot/plot_spacer()/secondary.hyp.g.moderator.failure.definition.plot/plot_spacer()/
  secondary.hyp.g.moderator.alternative.set.structure.plot/plot_spacer()/secondary.hyp.g.moderator.train.intent.plot&
  plot_annotation(title = "Moderator Analyses from Linear Multilevel Model for Muscle Hypertrophy (Standardized Mean Change)")&
  THEME

ggsave(secondary.hyp.g.moderator.doc,device = "pdf",filename="secondary.hyp.g.moderator.doc.pdf",width = 12,height = 136.5,limitsize = FALSE,
       path = "/Users/zacrobinson/Documents/For Publication PTF Meta/Plots")


#### rom load
secondary.hyp.rom.moderator.load=secondary.hyp.rom.best.model.prep%>%
  update(~. + avg.rir:load.set)%>%
  robust.rma.mv(cluster = secondary.data.hyp$study)

secondary.hyp.rom.moderator.load.emms=secondary.hyp.rom.moderator.load%>%
  emmprep(at=list(avg.rir=0:10,
                  load.set=c(30,60,90)))%>%
  emmeans(~avg.rir+load.set,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.hyp.rom.moderator.load)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

secondary.hyp.rom.moderator.load.contrast=emmprep(secondary.hyp.rom.moderator.load,
                                                at=list(avg.rir=secondary.hyp.mean.rir,
                                                        load.set=secondary.hyp.mean.load.set))%>%
  emmeans(~avg.rir+load.set,
          weights = "prop")%>%
  contrast("consec",interaction=TRUE)%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

secondary.hyp.rom.moderator.load.contrast.label=paste0('Interaction Contrast (Load) = ',
                                                     bquote(.(round(secondary.hyp.rom.moderator.load.contrast[[3]],5))),
                                                     " [95% CI: ",bquote(.(round(secondary.hyp.rom.moderator.load.contrast[[6]],5))),
                                                     ", ",bquote(.(round(secondary.hyp.rom.moderator.load.contrast[[7]],5))),"]")

secondary.hyp.rom.moderator.load.plot=ggplot(data = secondary.data.hyp%>%
                                             mutate(load.set=ifelse(load.set<45,"30% of 1RM",
                                                                    ifelse(load.set>45&load.set<75,"60% of 1RM","90% of 1RM"))),
                                           aes(x=avg.rir,
                                               y=yi.rom.exp))+THEME+theme(plot.caption = element_text(face = "italic"))+
  facet_wrap(~load.set)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Load (% of 1RM) x <10 Estimated RIR",
       caption = secondary.hyp.rom.moderator.load.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = secondary.hyp.rom.moderator.load.emms%>%
                mutate(load.set=ifelse(load.set==30,"30% of 1RM",
                                       ifelse(load.set==60,"60% of 1RM","90% of 1RM"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = secondary.hyp.rom.moderator.load.emms%>%
                mutate(load.set=ifelse(load.set==30,"30% of 1RM",
                                       ifelse(load.set==60,"60% of 1RM","90% of 1RM"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = secondary.hyp.rom.moderator.load.emms%>%
              mutate(load.set=ifelse(load.set==30,"30% of 1RM",
                                     ifelse(load.set==60,"60% of 1RM","90% of 1RM"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### rom volume equating method
secondary.hyp.rom.moderator.volume.method=secondary.hyp.rom.best.model.prep%>%
  update(~. + avg.rir:set.rep.equated)%>%
  robust.rma.mv(cluster = secondary.data.hyp$study)

secondary.hyp.rom.moderator.volume.method.emms=secondary.hyp.rom.moderator.volume.method%>%emmprep(at=list(avg.rir=0:10))%>%
  emmeans(~avg.rir+set.rep.equated,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.hyp.rom.moderator.volume.method)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

secondary.hyp.rom.moderator.volume.method.contrast=emmprep(secondary.hyp.rom.moderator.volume.method,
                                                         at=list(avg.rir=secondary.hyp.mean.rir))%>%
  emmeans(~avg.rir+set.rep.equated,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

secondary.hyp.rom.moderator.volume.method.contrast.label=paste0('Interaction Contrast (Set - Rep) = ',
                                                              bquote(.(round(secondary.hyp.rom.moderator.volume.method.contrast[[1,3]],5))),
                                                              " [95% CI: ",bquote(.(round(secondary.hyp.rom.moderator.volume.method.contrast[[1,6]],5))),
                                                              ", ",bquote(.(round(secondary.hyp.rom.moderator.volume.method.contrast[[1,7]],5))),"]\n",
                                                              'Interaction Contrast (Set - Both) = ',
                                                              bquote(.(round(secondary.hyp.rom.moderator.volume.method.contrast[[2,3]],5))),
                                                              " [95% CI: ",bquote(.(round(secondary.hyp.rom.moderator.volume.method.contrast[[2,6]],5))),
                                                              ", ",bquote(.(round(secondary.hyp.rom.moderator.volume.method.contrast[[2,7]],5))),"]\n",
                                                              'Interaction Contrast (Rep - Both) = ',
                                                              bquote(.(round(secondary.hyp.rom.moderator.volume.method.contrast[[3,3]],5))),
                                                              " [95% CI: ",bquote(.(round(secondary.hyp.rom.moderator.volume.method.contrast[[3,6]],5))),
                                                              ", ",bquote(.(round(secondary.hyp.rom.moderator.volume.method.contrast[[3,7]],5))),"]")


secondary.hyp.rom.moderator.volume.method.plot=ggplot(data = secondary.data.hyp%>%
                                                      mutate(set.rep.equated=ifelse(set.rep.equated=="set","Set Equated",
                                                                                    ifelse(set.rep.equated=="rep","Repetition Equated","Both Equated"))),
                                                    aes(x=avg.rir,
                                                        y=yi.rom.exp))+THEME+theme(legend.position = "right",
                                                                                   plot.caption = element_text(face = "italic"))+
  facet_wrap(~set.rep.equated)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Method of Volume Equating x <10 Estimated RIR",
       caption = secondary.hyp.rom.moderator.volume.method.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = secondary.hyp.rom.moderator.volume.method.emms%>%
                mutate(set.rep.equated=ifelse(set.rep.equated=="set","Set Equated",
                                              ifelse(set.rep.equated=="rep","Repetition Equated","Both Equated"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = secondary.hyp.rom.moderator.volume.method.emms%>%
                mutate(set.rep.equated=ifelse(set.rep.equated=="set","Set Equated",
                                              ifelse(set.rep.equated=="rep","Repetition Equated","Both Equated"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = secondary.hyp.rom.moderator.volume.method.emms%>%
              mutate(set.rep.equated=ifelse(set.rep.equated=="set","Set Equated",
                                            ifelse(set.rep.equated=="rep","Repetition Equated","Both Equated"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### rom training status
secondary.hyp.rom.moderator.train.status=secondary.hyp.rom.best.model.prep%>%
  update(~. + avg.rir:train.status)%>%
  robust.rma.mv(cluster = secondary.data.hyp$study)

secondary.hyp.rom.moderator.train.status.emms=secondary.hyp.rom.moderator.train.status%>%emmprep(at=list(avg.rir=0:10))%>%
  emmeans(~avg.rir+train.status,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.hyp.rom.moderator.train.status)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

secondary.hyp.rom.moderator.train.status.contrast=emmprep(secondary.hyp.rom.moderator.train.status,
                                                        at=list(avg.rir=secondary.hyp.mean.rir))%>%
  emmeans(~avg.rir+train.status,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

secondary.hyp.rom.moderator.train.status.contrast.label=paste0('Interaction Contrast (Trained - Untrained) = ',
                                                             bquote(.(round(secondary.hyp.rom.moderator.train.status.contrast[[1,3]],5))),
                                                             " [95% CI: ",bquote(.(round(secondary.hyp.rom.moderator.train.status.contrast[[1,6]],5))),
                                                             ", ",bquote(.(round(secondary.hyp.rom.moderator.train.status.contrast[[1,7]],5))),"]")


secondary.hyp.rom.moderator.train.status.plot=ggplot(data = secondary.data.hyp%>%
                                                     mutate(train.status=ifelse(train.status=="trained","Trained","Untrained")),
                                                   aes(x=avg.rir,
                                                       y=yi.rom.exp))+THEME+theme(legend.position = "right",
                                                                                  plot.caption = element_text(face = "italic"))+
  facet_wrap(~train.status)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Training Status x <10 Estimated RIR",
       caption = secondary.hyp.rom.moderator.train.status.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = secondary.hyp.rom.moderator.train.status.emms%>%
                mutate(train.status=ifelse(train.status=="trained","Trained","Untrained")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = secondary.hyp.rom.moderator.train.status.emms%>%
                mutate(train.status=ifelse(train.status=="trained","Trained","Untrained")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = secondary.hyp.rom.moderator.train.status.emms%>%
              mutate(train.status=ifelse(train.status=="trained","Trained","Untrained")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### rom weeks
secondary.hyp.rom.moderator.weeks=secondary.hyp.rom.best.model.prep%>%
  update(~. + avg.rir:weeks)%>%
  robust.rma.mv(cluster = secondary.data.hyp$study)

secondary.hyp.rom.moderator.weeks.emms=secondary.hyp.rom.moderator.weeks%>%
  emmprep(at=list(avg.rir=0:10,
                  weeks=c(4,8,12)))%>%
  emmeans(~avg.rir+weeks,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.hyp.rom.moderator.weeks)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

secondary.hyp.rom.moderator.weeks.contrast=emmprep(secondary.hyp.rom.moderator.weeks,
                                                 at=list(avg.rir=secondary.hyp.mean.rir,
                                                         weeks=secondary.hyp.mean.weeks))%>%
  emmeans(~avg.rir+weeks,
          weights = "prop")%>%
  contrast("consec",interaction=TRUE)%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

secondary.hyp.rom.moderator.weeks.contrast.label=paste0('Interaction Contrast (Weeks) = ',
                                                      bquote(.(round(secondary.hyp.rom.moderator.weeks.contrast[[3]],5))),
                                                      " [95% CI: ",bquote(.(round(secondary.hyp.rom.moderator.weeks.contrast[[6]],5))),
                                                      ", ",bquote(.(round(secondary.hyp.rom.moderator.weeks.contrast[[7]],5))),"]")

secondary.hyp.rom.moderator.weeks.plot=ggplot(data = secondary.data.hyp%>%
                                              mutate(weeks=ifelse(weeks<6,"4 Weeks",
                                                                  ifelse(weeks>6&weeks<10,"8 Weeks","12 Weeks"))),
                                            aes(x=avg.rir,
                                                y=yi.rom.exp))+THEME+theme(plot.caption = element_text(face = "italic"))+
  facet_wrap(~fct_relevel(weeks,c("4 Weeks","8 Weeks","12 Weeks")))+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Intervention Duration (Weeks) x <10 Estimated RIR",
       caption = secondary.hyp.rom.moderator.weeks.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = secondary.hyp.rom.moderator.weeks.emms%>%
                mutate(weeks=ifelse(weeks==4,"4 Weeks",
                                    ifelse(weeks==8,"8 Weeks","12 Weeks"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = secondary.hyp.rom.moderator.weeks.emms%>%
                mutate(weeks=ifelse(weeks==4,"4 Weeks",
                                    ifelse(weeks==8,"8 Weeks","12 Weeks"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = secondary.hyp.rom.moderator.weeks.emms%>%
              mutate(weeks=ifelse(weeks==4,"4 Weeks",
                                  ifelse(weeks==8,"8 Weeks","12 Weeks"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### rom adjusted sets 
secondary.hyp.rom.moderator.adj.sets=secondary.hyp.rom.best.model.prep%>%
  update(~. + avg.rir*adj.sets.week)%>%
  robust.rma.mv(cluster = secondary.data.hyp$study)

secondary.hyp.rom.moderator.adj.sets.emms=secondary.hyp.rom.moderator.adj.sets%>%
  emmprep(at=list(avg.rir=0:10,
                  adj.sets.week=c(8,16,24)))%>%
  emmeans(~avg.rir+adj.sets.week,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.hyp.rom.moderator.adj.sets)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

secondary.hyp.rom.moderator.adj.sets.contrast=emmprep(secondary.hyp.rom.moderator.adj.sets,
                                                    at=list(avg.rir=secondary.hyp.mean.rir,
                                                            adj.sets.week=secondary.hyp.mean.adj.sets.week))%>%
  emmeans(~avg.rir+adj.sets.week,
          weights = "prop")%>%
  contrast("consec",interaction=TRUE)%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

secondary.hyp.rom.moderator.adj.sets.contrast.label=paste0('Interaction Contrast (Adjusted Direct Sets) = ',
                                                         bquote(.(round(secondary.hyp.rom.moderator.adj.sets.contrast[[3]],5))),
                                                         " [95% CI: ",bquote(.(round(secondary.hyp.rom.moderator.adj.sets.contrast[[6]],5))),
                                                         ", ",bquote(.(round(secondary.hyp.rom.moderator.adj.sets.contrast[[7]],5))),"]")

secondary.hyp.rom.moderator.adj.sets.plot=ggplot(data = secondary.data.hyp%>%
                                                 mutate(adj.sets.week=ifelse(adj.sets.week<12,"8 Sets",
                                                                             ifelse(adj.sets.week>12&adj.sets.week<20,"16 Sets","24 Sets"))),
                                               aes(x=avg.rir,
                                                   y=yi.rom.exp))+THEME+theme(plot.caption = element_text(face = "italic"))+
  facet_wrap(~fct_relevel(adj.sets.week,c("8 Sets","16 Sets","24 Sets")))+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Adjusted Direct Sets Per Week Per Muscle Group x <10 Estimated RIR",
       caption = secondary.hyp.rom.moderator.adj.sets.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = secondary.hyp.rom.moderator.adj.sets.emms%>%
                mutate(adj.sets.week=ifelse(adj.sets.week==8,"8 Sets",
                                            ifelse(adj.sets.week==16,"16 Sets","24 Sets"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = secondary.hyp.rom.moderator.adj.sets.emms%>%
                mutate(adj.sets.week=ifelse(adj.sets.week==8,"8 Sets",
                                            ifelse(adj.sets.week==16,"16 Sets","24 Sets"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = secondary.hyp.rom.moderator.adj.sets.emms%>%
              mutate(adj.sets.week=ifelse(adj.sets.week==8,"8 Sets",
                                          ifelse(adj.sets.week==16,"16 Sets","24 Sets"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### rom frequency
secondary.hyp.rom.moderator.frequency.per.muscle=secondary.hyp.rom.best.model.prep%>%
  update(~. + avg.rir*frequency.per.muscle)%>%
  robust.rma.mv(cluster = secondary.data.hyp$study)

secondary.hyp.rom.moderator.frequency.per.muscle.emms=secondary.hyp.rom.moderator.frequency.per.muscle%>%
  emmprep(at=list(avg.rir=0:10,
                  frequency.per.muscle=c(1,2,3)))%>%
  emmeans(~avg.rir+frequency.per.muscle,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.hyp.rom.moderator.frequency.per.muscle)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

secondary.hyp.rom.moderator.frequency.per.muscle.contrast=emmprep(secondary.hyp.rom.moderator.frequency.per.muscle,
                                                                at=list(avg.rir=secondary.hyp.mean.rir,
                                                                        frequency.per.muscle=secondary.hyp.mean.frequency.per.muscle))%>%
  emmeans(~avg.rir+frequency.per.muscle,
          weights = "prop")%>%
  contrast("consec",interaction=TRUE)%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

secondary.hyp.rom.moderator.frequency.per.muscle.contrast.label=paste0('Interaction Contrast (Direct Frequency) = ',
                                                                     bquote(.(round(secondary.hyp.rom.moderator.frequency.per.muscle.contrast[[3]],5))),
                                                                     " [95% CI: ",bquote(.(round(secondary.hyp.rom.moderator.frequency.per.muscle.contrast[[6]],5))),
                                                                     ", ",bquote(.(round(secondary.hyp.rom.moderator.frequency.per.muscle.contrast[[7]],5))),"]")

secondary.hyp.rom.moderator.frequency.per.muscle.plot=ggplot(data = secondary.data.hyp%>%
                                                             mutate(frequency.per.muscle=ifelse(frequency.per.muscle<1.5,"1x Per Week",
                                                                                                ifelse(frequency.per.muscle>1.5&frequency.per.muscle<2.5,"2x Per Week","3x Per Week"))),
                                                           aes(x=avg.rir,
                                                               y=yi.rom.exp))+THEME+theme(plot.caption = element_text(face = "italic"))+
  facet_wrap(~frequency.per.muscle)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Direct Frequency Per Week Per Muscle Group x <10 Estimated RIR",
       caption = secondary.hyp.rom.moderator.frequency.per.muscle.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = secondary.hyp.rom.moderator.frequency.per.muscle.emms%>%
                mutate(frequency.per.muscle=ifelse(frequency.per.muscle==1,"1x Per Week",
                                                   ifelse(frequency.per.muscle==2,"2x Per Week","3x Per Week"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = secondary.hyp.rom.moderator.frequency.per.muscle.emms%>%
                mutate(frequency.per.muscle=ifelse(frequency.per.muscle==1,"1x Per Week",
                                                   ifelse(frequency.per.muscle==2,"2x Per Week","3x Per Week"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = secondary.hyp.rom.moderator.frequency.per.muscle.emms%>%
              mutate(frequency.per.muscle=ifelse(frequency.per.muscle==1,"1x Per Week",
                                                 ifelse(frequency.per.muscle==2,"2x Per Week","3x Per Week"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### rom age
secondary.hyp.rom.moderator.age=secondary.hyp.rom.best.model.prep%>%
  update(~. + avg.rir*age)%>%
  robust.rma.mv(cluster = secondary.data.hyp$study)

secondary.hyp.rom.moderator.age.emms=secondary.hyp.rom.moderator.age%>%
  emmprep(at=list(avg.rir=0:10,
                  age=c(20,40,60)))%>%
  emmeans(~avg.rir+age,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.hyp.rom.moderator.age)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

secondary.hyp.rom.moderator.age.contrast=emmprep(secondary.hyp.rom.moderator.age,
                                               at=list(avg.rir=secondary.hyp.mean.rir,
                                                       age=secondary.hyp.mean.age))%>%
  emmeans(~avg.rir+age,
          weights = "prop")%>%
  contrast("consec",interaction=TRUE)%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

secondary.hyp.rom.moderator.age.contrast.label=paste0('Interaction Contrast (Age) = ',
                                                    bquote(.(round(secondary.hyp.rom.moderator.age.contrast[[3]],5))),
                                                    " [95% CI: ",bquote(.(round(secondary.hyp.rom.moderator.age.contrast[[6]],5))),
                                                    ", ",bquote(.(round(secondary.hyp.rom.moderator.age.contrast[[7]],5))),"]")

secondary.hyp.rom.moderator.age.plot=ggplot(data = secondary.data.hyp%>%
                                            mutate(age=ifelse(age<30,"20 Years Old",
                                                              ifelse(age>30&age<50,"40 Years Old","60 Years Old"))),
                                          aes(x=avg.rir,
                                              y=yi.rom.exp))+THEME+theme(plot.caption = element_text(face = "italic"))+
  facet_wrap(~age)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Age x <10 Estimated RIR",
       caption = secondary.hyp.rom.moderator.age.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = secondary.hyp.rom.moderator.age.emms%>%
                mutate(age=ifelse(age==20,"20 Years Old",
                                  ifelse(age==40,"40 Years Old","60 Years Old"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = secondary.hyp.rom.moderator.age.emms%>%
                mutate(age=ifelse(age==20,"20 Years Old",
                                  ifelse(age==40,"40 Years Old","60 Years Old"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = secondary.hyp.rom.moderator.age.emms%>%
              mutate(age=ifelse(age==20,"20 Years Old",
                                ifelse(age==40,"40 Years Old","60 Years Old"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### rom sex
secondary.hyp.rom.moderator.sex.percent.male=secondary.hyp.rom.best.model.prep%>%
  update(~. + avg.rir*sex.percent.male)%>%
  robust.rma.mv(cluster = secondary.data.hyp$study)

secondary.hyp.rom.moderator.sex.percent.male.emms=secondary.hyp.rom.moderator.sex.percent.male%>%
  emmprep(at=list(avg.rir=0:10,
                  sex.percent.male=c(0,100)))%>%
  emmeans(~avg.rir+sex.percent.male,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.hyp.rom.moderator.sex.percent.male)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

secondary.hyp.rom.moderator.sex.percent.male.contrast=emmprep(secondary.hyp.rom.moderator.sex.percent.male,
                                                            at=list(avg.rir=secondary.hyp.mean.rir,
                                                                    sex.percent.male=secondary.hyp.mean.sex))%>%
  emmeans(~avg.rir+sex.percent.male,
          weights = "prop")%>%
  contrast("consec",interaction=TRUE)%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

secondary.hyp.rom.moderator.sex.percent.male.contrast.label=paste0('Interaction Contrast (Sex) = ',
                                                                 bquote(.(round(secondary.hyp.rom.moderator.sex.percent.male.contrast[[3]],5))),
                                                                 " [95% CI: ",bquote(.(round(secondary.hyp.rom.moderator.sex.percent.male.contrast[[6]],5))),
                                                                 ", ",bquote(.(round(secondary.hyp.rom.moderator.sex.percent.male.contrast[[7]],5))),"]")

secondary.hyp.rom.moderator.sex.percent.male.plot=ggplot(data = secondary.data.hyp%>%
                                                         mutate(sex.percent.male=ifelse(sex.percent.male<50,"Female","Male")),
                                                       aes(x=avg.rir,
                                                           y=yi.rom.exp))+THEME+theme(plot.caption = element_text(face = "italic"))+
  facet_wrap(~sex.percent.male)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Sex (% Male) x <10 Estimated RIR",
       caption = secondary.hyp.rom.moderator.sex.percent.male.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = secondary.hyp.rom.moderator.sex.percent.male.emms%>%
                mutate(sex.percent.male=ifelse(sex.percent.male==0,"Female","Male")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = secondary.hyp.rom.moderator.sex.percent.male.emms%>%
                mutate(sex.percent.male=ifelse(sex.percent.male==0,"Female","Male")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = secondary.hyp.rom.moderator.sex.percent.male.emms%>%
              mutate(sex.percent.male=ifelse(sex.percent.male==0,"Female","Male")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### rom study design
secondary.hyp.rom.moderator.within.between.design=secondary.hyp.rom.best.model.prep%>%
  update(~. + avg.rir*within.between.design)%>%
  robust.rma.mv(cluster = secondary.data.hyp$study)

secondary.hyp.rom.moderator.within.between.design.emms=secondary.hyp.rom.moderator.within.between.design%>%emmprep(at=list(avg.rir=0:10))%>%
  emmeans(~avg.rir+within.between.design,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.hyp.rom.moderator.within.between.design)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

secondary.hyp.rom.moderator.within.between.design.contrast=emmprep(secondary.hyp.rom.moderator.within.between.design,
                                                                 at=list(avg.rir=secondary.hyp.mean.rir))%>%
  emmeans(~avg.rir+within.between.design,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

secondary.hyp.rom.moderator.within.between.design.contrast.label=paste0('Interaction Contrast (Between - Within) = ',
                                                                      bquote(.(round(secondary.hyp.rom.moderator.within.between.design.contrast[[1,3]],5))),
                                                                      " [95% CI: ",bquote(.(round(secondary.hyp.rom.moderator.within.between.design.contrast[[1,6]],5))),
                                                                      ", ",bquote(.(round(secondary.hyp.rom.moderator.within.between.design.contrast[[1,7]],5))),"]")


secondary.hyp.rom.moderator.within.between.design.plot=ggplot(data = secondary.data.hyp%>%
                                                              mutate(within.between.design=ifelse(within.between.design=="Between","Between Subjects","Within Subjects")),
                                                            aes(x=avg.rir,
                                                                y=yi.rom.exp))+THEME+theme(legend.position = "right",
                                                                                           plot.caption = element_text(face = "italic"))+
  facet_wrap(~within.between.design)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Study Design x <10 Estimated RIR",
       caption = secondary.hyp.rom.moderator.within.between.design.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = secondary.hyp.rom.moderator.within.between.design.emms%>%
                mutate(within.between.design=ifelse(within.between.design=="Between","Between Subjects","Within Subjects")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = secondary.hyp.rom.moderator.within.between.design.emms%>%
                mutate(within.between.design=ifelse(within.between.design=="Between","Between Subjects","Within Subjects")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = secondary.hyp.rom.moderator.within.between.design.emms%>%
              mutate(within.between.design=ifelse(within.between.design=="Between","Between Subjects","Within Subjects")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### rom upper / lower outcome
secondary.hyp.rom.moderator.upper.lower.other=secondary.hyp.rom.best.model.prep%>%
  update(~. + avg.rir*upper.lower.other)%>%
  robust.rma.mv(cluster = secondary.data.hyp$study)

secondary.hyp.rom.moderator.upper.lower.other.emms=secondary.hyp.rom.moderator.upper.lower.other%>%emmprep(at=list(avg.rir=0:10))%>%
  emmeans(~avg.rir+upper.lower.other,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.hyp.rom.moderator.upper.lower.other)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

secondary.hyp.rom.moderator.upper.lower.other.contrast=emmprep(secondary.hyp.rom.moderator.upper.lower.other,
                                                             at=list(avg.rir=secondary.hyp.mean.rir))%>%
  emmeans(~avg.rir+upper.lower.other,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

secondary.hyp.rom.moderator.upper.lower.other.contrast.label=paste0('Interaction Contrast (Lower - Upper) = ',
                                                                  bquote(.(round(secondary.hyp.rom.moderator.upper.lower.other.contrast[[1,3]],5))),
                                                                  " [95% CI: ",bquote(.(round(secondary.hyp.rom.moderator.upper.lower.other.contrast[[1,6]],5))),
                                                                  ", ",bquote(.(round(secondary.hyp.rom.moderator.upper.lower.other.contrast[[1,7]],5))),"]")


secondary.hyp.rom.moderator.upper.lower.other.plot=ggplot(data = secondary.data.hyp%>%
                                                          mutate(upper.lower.other=ifelse(upper.lower.other=="Lower","Lower Body Outcomes","Upper Body Outcomes")),
                                                        aes(x=avg.rir,
                                                            y=yi.rom.exp))+THEME+theme(legend.position = "right",
                                                                                       plot.caption = element_text(face = "italic"))+
  facet_wrap(~upper.lower.other)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Upper/Lower Body Outcome x <10 Estimated RIR",
       caption = secondary.hyp.rom.moderator.upper.lower.other.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = secondary.hyp.rom.moderator.upper.lower.other.emms%>%
                mutate(upper.lower.other=ifelse(upper.lower.other=="Lower","Lower Body Outcomes","Upper Body Outcomes")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = secondary.hyp.rom.moderator.upper.lower.other.emms%>%
                mutate(upper.lower.other=ifelse(upper.lower.other=="Lower","Lower Body Outcomes","Upper Body Outcomes")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = secondary.hyp.rom.moderator.upper.lower.other.emms%>%
              mutate(upper.lower.other=ifelse(upper.lower.other=="Lower","Lower Body Outcomes","Upper Body Outcomes")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### rom exercise selection
secondary.hyp.rom.moderator.train.exercise=secondary.hyp.rom.best.model.prep%>%
  update(~. + avg.rir*train.exercise)%>%
  robust.rma.mv(cluster = secondary.data.hyp$study)

secondary.hyp.rom.moderator.train.exercise.emms=secondary.hyp.rom.moderator.train.exercise%>%emmprep(at=list(avg.rir=0:10))%>%
  emmeans(~avg.rir+train.exercise,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.hyp.rom.moderator.train.exercise)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

secondary.hyp.rom.moderator.train.exercise.contrast=emmprep(secondary.hyp.rom.moderator.train.exercise,
                                                          at=list(avg.rir=secondary.hyp.mean.rir))%>%
  emmeans(~avg.rir+train.exercise,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

secondary.hyp.rom.moderator.train.exercise.contrast.label=paste0('Interaction Contrast (Both - Multi-joint) = ',
                                                               bquote(.(round(secondary.hyp.rom.moderator.train.exercise.contrast[[1,3]],5))),
                                                               " [95% CI: ",bquote(.(round(secondary.hyp.rom.moderator.train.exercise.contrast[[1,6]],5))),
                                                               ", ",bquote(.(round(secondary.hyp.rom.moderator.train.exercise.contrast[[1,7]],5))),"]\n",
                                                               'Interaction Contrast (Both - Single-joint) = ',
                                                               bquote(.(round(secondary.hyp.rom.moderator.train.exercise.contrast[[2,3]],5))),
                                                               " [95% CI: ",bquote(.(round(secondary.hyp.rom.moderator.train.exercise.contrast[[2,6]],5))),
                                                               ", ",bquote(.(round(secondary.hyp.rom.moderator.train.exercise.contrast[[2,7]],5))),"]\n",
                                                               'Interaction Contrast (Multi-joint - Single-joint) = ',
                                                               bquote(.(round(secondary.hyp.rom.moderator.train.exercise.contrast[[3,3]],5))),
                                                               " [95% CI: ",bquote(.(round(secondary.hyp.rom.moderator.train.exercise.contrast[[3,6]],5))),
                                                               ", ",bquote(.(round(secondary.hyp.rom.moderator.train.exercise.contrast[[3,7]],5))),"]")


secondary.hyp.rom.moderator.train.exercise.plot=ggplot(data = secondary.data.hyp%>%
                                                       mutate(train.exercise=ifelse(train.exercise=="Both","Both",
                                                                                    ifelse(train.exercise=="Multi-joint","Multi-joint Exercises","Single-joint Exercises"))),
                                                     aes(x=avg.rir,
                                                         y=yi.rom.exp))+THEME+theme(legend.position = "right",
                                                                                    plot.caption = element_text(face = "italic"))+
  facet_wrap(~train.exercise)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Exercise Selection x <10 Estimated RIR",
       caption = secondary.hyp.rom.moderator.train.exercise.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = secondary.hyp.rom.moderator.train.exercise.emms%>%
                mutate(train.exercise=ifelse(train.exercise=="Both","Both",
                                             ifelse(train.exercise=="Multi-joint","Multi-joint Exercises","Single-joint Exercises"))),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = secondary.hyp.rom.moderator.train.exercise.emms%>%
                mutate(train.exercise=ifelse(train.exercise=="Both","Both",
                                             ifelse(train.exercise=="Multi-joint","Multi-joint Exercises","Single-joint Exercises"))),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = secondary.hyp.rom.moderator.train.exercise.emms%>%
              mutate(train.exercise=ifelse(train.exercise=="Both","Both",
                                           ifelse(train.exercise=="Multi-joint","Multi-joint Exercises","Single-joint Exercises"))),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### rom concurrent training
secondary.hyp.rom.moderator.formal.cardio.intervention=secondary.hyp.rom.best.model.prep%>%
  update(~. + avg.rir*formal.cardio.intervention)%>%
  robust.rma.mv(cluster = secondary.data.hyp$study)

secondary.hyp.rom.moderator.formal.cardio.intervention.emms=secondary.hyp.rom.moderator.formal.cardio.intervention%>%emmprep(at=list(avg.rir=0:10))%>%
  emmeans(~avg.rir+formal.cardio.intervention,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.hyp.rom.moderator.formal.cardio.intervention)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

secondary.hyp.rom.moderator.formal.cardio.intervention.contrast=emmprep(secondary.hyp.rom.moderator.formal.cardio.intervention,
                                                                      at=list(avg.rir=secondary.hyp.mean.rir))%>%
  emmeans(~avg.rir+formal.cardio.intervention,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

secondary.hyp.rom.moderator.formal.cardio.intervention.contrast.label=paste0('Interaction Contrast (RT Only - Concurrent) = ',
                                                                           bquote(.(round(secondary.hyp.rom.moderator.formal.cardio.intervention.contrast[[1,3]],5))),
                                                                           " [95% CI: ",bquote(.(round(secondary.hyp.rom.moderator.formal.cardio.intervention.contrast[[1,6]],5))),
                                                                           ", ",bquote(.(round(secondary.hyp.rom.moderator.formal.cardio.intervention.contrast[[1,7]],5))),"]")


secondary.hyp.rom.moderator.formal.cardio.intervention.plot=ggplot(data = secondary.data.hyp%>%
                                                                   mutate(formal.cardio.intervention=ifelse(formal.cardio.intervention=="No","Only Resistance Training","Concurrent Training")),
                                                                 aes(x=avg.rir,
                                                                     y=yi.rom.exp))+THEME+theme(legend.position = "right",
                                                                                                plot.caption = element_text(face = "italic"))+
  facet_wrap(~formal.cardio.intervention)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Concurrent Training x <10 Estimated RIR",
       caption = secondary.hyp.rom.moderator.formal.cardio.intervention.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = secondary.hyp.rom.moderator.formal.cardio.intervention.emms%>%
                mutate(formal.cardio.intervention=ifelse(formal.cardio.intervention=="No","Only Resistance Training","Concurrent Training")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = secondary.hyp.rom.moderator.formal.cardio.intervention.emms%>%
                mutate(formal.cardio.intervention=ifelse(formal.cardio.intervention=="No","Only Resistance Training","Concurrent Training")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = secondary.hyp.rom.moderator.formal.cardio.intervention.emms%>%
              mutate(formal.cardio.intervention=ifelse(formal.cardio.intervention=="No","Only Resistance Training","Concurrent Training")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### rom progressive overload
secondary.hyp.rom.moderator.progression=secondary.hyp.rom.best.model.prep%>%
  update(~. + avg.rir*progression)%>%
  robust.rma.mv(cluster = secondary.data.hyp$study)

secondary.hyp.rom.moderator.progression.emms=secondary.hyp.rom.moderator.progression%>%emmprep(at=list(avg.rir=0:10))%>%
  emmeans(~avg.rir+progression,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.hyp.rom.moderator.progression)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

secondary.hyp.rom.moderator.progression.contrast=emmprep(secondary.hyp.rom.moderator.progression,
                                                       at=list(avg.rir=secondary.hyp.mean.rir))%>%
  emmeans(~avg.rir+progression,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

secondary.hyp.rom.moderator.progression.contrast.label=paste0('Interaction Contrast (No Overload - Progressive Overload) = ',
                                                            bquote(.(round(secondary.hyp.rom.moderator.progression.contrast[[1,3]],5))),
                                                            " [95% CI: ",bquote(.(round(secondary.hyp.rom.moderator.progression.contrast[[1,6]],5))),
                                                            ", ",bquote(.(round(secondary.hyp.rom.moderator.progression.contrast[[1,7]],5))),"]")


secondary.hyp.rom.moderator.progression.plot=ggplot(data = secondary.data.hyp%>%
                                                    mutate(progression=ifelse(progression=="No","No Progressive Overload","Progressive Overload")),
                                                  aes(x=avg.rir,
                                                      y=yi.rom.exp))+THEME+theme(legend.position = "right",
                                                                                 plot.caption = element_text(face = "italic"))+
  facet_wrap(~progression)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Progressive Overload x <10 Estimated RIR",
       caption = secondary.hyp.rom.moderator.progression.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = secondary.hyp.rom.moderator.progression.emms%>%
                mutate(progression=ifelse(progression=="No","No Progressive Overload","Progressive Overload")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = secondary.hyp.rom.moderator.progression.emms%>%
                mutate(progression=ifelse(progression=="No","No Progressive Overload","Progressive Overload")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = secondary.hyp.rom.moderator.progression.emms%>%
              mutate(progression=ifelse(progression=="No","No Progressive Overload","Progressive Overload")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### rom failure definition
secondary.hyp.rom.moderator.failure.definition=secondary.hyp.rom.best.model.prep%>%
  update(~. + avg.rir*failure.definition)%>%
  robust.rma.mv(cluster = secondary.data.hyp$study)

secondary.hyp.rom.moderator.failure.definition.emms=secondary.hyp.rom.moderator.failure.definition%>%emmprep(at=list(avg.rir=0:10))%>%
  emmeans(~avg.rir+failure.definition,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.hyp.rom.moderator.failure.definition)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

secondary.hyp.rom.moderator.failure.definition.contrast=emmprep(secondary.hyp.rom.moderator.failure.definition,
                                                              at=list(avg.rir=secondary.hyp.mean.rir))%>%
  emmeans(~avg.rir+failure.definition,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

secondary.hyp.rom.moderator.failure.definition.contrast.label=paste0('Interaction Contrast (No Definition - Definition Provided) = ',
                                                                   bquote(.(round(secondary.hyp.rom.moderator.failure.definition.contrast[[1,3]],5))),
                                                                   " [95% CI: ",bquote(.(round(secondary.hyp.rom.moderator.failure.definition.contrast[[1,6]],5))),
                                                                   ", ",bquote(.(round(secondary.hyp.rom.moderator.failure.definition.contrast[[1,7]],5))),"]")


secondary.hyp.rom.moderator.failure.definition.plot=ggplot(data = secondary.data.hyp%>%
                                                           mutate(failure.definition=ifelse(failure.definition=="No","No Failure Definition","Failure Definition Provided")),
                                                         aes(x=avg.rir,
                                                             y=yi.rom.exp))+THEME+theme(legend.position = "right",
                                                                                        plot.caption = element_text(face = "italic"))+
  facet_wrap(~failure.definition)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Failure Definition x <10 Estimated RIR",
       caption = secondary.hyp.rom.moderator.failure.definition.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = secondary.hyp.rom.moderator.failure.definition.emms%>%
                mutate(failure.definition=ifelse(failure.definition=="No","No Failure Definition","Failure Definition Provided")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = secondary.hyp.rom.moderator.failure.definition.emms%>%
                mutate(failure.definition=ifelse(failure.definition=="No","No Failure Definition","Failure Definition Provided")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = secondary.hyp.rom.moderator.failure.definition.emms%>%
              mutate(failure.definition=ifelse(failure.definition=="No","No Failure Definition","Failure Definition Provided")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### rom set structure
secondary.hyp.rom.moderator.alternative.set.structure=secondary.hyp.rom.best.model.prep%>%
  update(~. + avg.rir*alternative.set.structure)%>%
  robust.rma.mv(cluster = secondary.data.hyp$study)

secondary.hyp.rom.moderator.alternative.set.structure.emms=secondary.hyp.rom.moderator.alternative.set.structure%>%emmprep(at=list(avg.rir=0:10))%>%
  emmeans(~avg.rir+alternative.set.structure,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.hyp.rom.moderator.alternative.set.structure)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

secondary.hyp.rom.moderator.alternative.set.structure.contrast=emmprep(secondary.hyp.rom.moderator.alternative.set.structure,
                                                                     at=list(avg.rir=secondary.hyp.mean.rir))%>%
  emmeans(~avg.rir+alternative.set.structure,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

secondary.hyp.rom.moderator.alternative.set.structure.contrast.label=paste0('Interaction Contrast (Traditional - Alternative) = ',
                                                                          bquote(.(round(secondary.hyp.rom.moderator.alternative.set.structure.contrast[[1,3]],5))),
                                                                          " [95% CI: ",bquote(.(round(secondary.hyp.rom.moderator.alternative.set.structure.contrast[[1,6]],5))),
                                                                          ", ",bquote(.(round(secondary.hyp.rom.moderator.alternative.set.structure.contrast[[1,7]],5))),"]")


secondary.hyp.rom.moderator.alternative.set.structure.plot=ggplot(data = secondary.data.hyp%>%
                                                                  mutate(alternative.set.structure=ifelse(alternative.set.structure=="No","Traditional Set Structure","Alternative Set Structure")),
                                                                aes(x=avg.rir,
                                                                    y=yi.rom.exp))+THEME+theme(legend.position = "right",
                                                                                               plot.caption = element_text(face = "italic"))+
  facet_wrap(~alternative.set.structure)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Set Structure x <10 Estimated RIR",
       caption = secondary.hyp.rom.moderator.alternative.set.structure.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = secondary.hyp.rom.moderator.alternative.set.structure.emms%>%
                mutate(alternative.set.structure=ifelse(alternative.set.structure=="No","Traditional Set Structure","Alternative Set Structure")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = secondary.hyp.rom.moderator.alternative.set.structure.emms%>%
                mutate(alternative.set.structure=ifelse(alternative.set.structure=="No","Traditional Set Structure","Alternative Set Structure")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = secondary.hyp.rom.moderator.alternative.set.structure.emms%>%
              mutate(alternative.set.structure=ifelse(alternative.set.structure=="No","Traditional Set Structure","Alternative Set Structure")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### rom intended velocity
secondary.hyp.rom.moderator.train.intent=secondary.hyp.rom.best.model.prep%>%
  update(~. + avg.rir*train.intent)%>%
  robust.rma.mv(cluster = secondary.data.hyp$study)

secondary.hyp.rom.moderator.train.intent.emms=secondary.hyp.rom.moderator.train.intent%>%emmprep(at=list(avg.rir=0:10))%>%
  emmeans(~avg.rir+train.intent,
          weights = "prop")%>%
  pred_interval_esmeans(model=secondary.hyp.rom.moderator.train.intent)%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7,8,9),expo)

secondary.hyp.rom.moderator.train.intent.contrast=emmprep(secondary.hyp.rom.moderator.train.intent,
                                                        at=list(avg.rir=secondary.hyp.mean.rir))%>%
  emmeans(~avg.rir+train.intent,
          weights = "prop")%>%
  contrast(interaction=c("consec","pairwise"))%>%
  confint()%>%
  as.data.frame()%>%
  mutate_at(c(3,4,6,7),expo)

secondary.hyp.rom.moderator.train.intent.contrast.label=paste0('Interaction Contrast (Submaximal - Maximal) = ',
                                                             bquote(.(round(secondary.hyp.rom.moderator.train.intent.contrast[[1,3]],5))),
                                                             " [95% CI: ",bquote(.(round(secondary.hyp.rom.moderator.train.intent.contrast[[1,6]],5))),
                                                             ", ",bquote(.(round(secondary.hyp.rom.moderator.train.intent.contrast[[1,7]],5))),"]")


secondary.hyp.rom.moderator.train.intent.plot=ggplot(data = secondary.data.hyp%>%
                                                     mutate(train.intent=ifelse(train.intent=="No","Submaximal Intended Concentric Velocity","Maximal Intended Concentric Velocity")),
                                                   aes(x=avg.rir,
                                                       y=yi.rom.exp))+THEME+theme(legend.position = "right",
                                                                                  plot.caption = element_text(face = "italic"))+
  facet_wrap(~train.intent)+
  labs(x="Estimated Proximity to Failure (RIR)",
       y="Exponentiated Response Ratio (% Change)",
       title = "Intended Concentric Velocity x <10 Estimated RIR",
       caption = secondary.hyp.rom.moderator.train.intent.contrast.label)+
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  geom_point(aes(size=weights),
             color="#0000a7",
             alpha=0.75)+
  geom_ribbon(data = secondary.hyp.rom.moderator.train.intent.emms%>%
                mutate(train.intent=ifelse(train.intent=="No","Submaximal Intended Concentric Velocity","Maximal Intended Concentric Velocity")),
              aes(y=emmean,ymin=lower.CL,ymax=upper.CL),alpha=0.25,fill="#0000a7")+
  geom_ribbon(data = secondary.hyp.rom.moderator.train.intent.emms%>%
                mutate(train.intent=ifelse(train.intent=="No","Submaximal Intended Concentric Velocity","Maximal Intended Concentric Velocity")),
              aes(y=emmean,ymin=lower.PI,ymax=upper.PI),fill="#0000a7",alpha=0.15)+
  geom_line(data = secondary.hyp.rom.moderator.train.intent.emms%>%
              mutate(train.intent=ifelse(train.intent=="No","Submaximal Intended Concentric Velocity","Maximal Intended Concentric Velocity")),
            aes(x=avg.rir,
                y=emmean),
            linewidth=1,
            color="#0000a7")+
  scale_size_continuous(guide="none")

#### combine rom plots
secondary.hyp.rom.moderator.doc=secondary.hyp.rom.moderator.load.plot/plot_spacer()/secondary.hyp.rom.moderator.volume.method.plot/plot_spacer()/
  secondary.hyp.rom.moderator.train.status.plot/plot_spacer()/secondary.hyp.rom.moderator.weeks.plot/plot_spacer()/
  secondary.hyp.rom.moderator.adj.sets.plot/plot_spacer()/secondary.hyp.rom.moderator.frequency.per.muscle.plot/plot_spacer()/
  secondary.hyp.rom.moderator.age.plot/plot_spacer()/secondary.hyp.rom.moderator.sex.percent.male.plot/plot_spacer()/
  secondary.hyp.rom.moderator.within.between.design.plot/plot_spacer()/secondary.hyp.rom.moderator.upper.lower.other.plot/plot_spacer()/
  secondary.hyp.rom.moderator.train.exercise.plot/plot_spacer()/secondary.hyp.rom.moderator.formal.cardio.intervention.plot/plot_spacer()/
  secondary.hyp.rom.moderator.progression.plot/plot_spacer()/secondary.hyp.rom.moderator.failure.definition.plot/plot_spacer()/
  secondary.hyp.rom.moderator.alternative.set.structure.plot/plot_spacer()/secondary.hyp.rom.moderator.train.intent.plot&
  plot_annotation(title = "Moderator Analyses from Linear Multilevel Model for Muscle Hypertrophy (Response Ratio)")&
  THEME

ggsave(secondary.hyp.rom.moderator.doc,device = "pdf",filename="secondary.hyp.rom.moderator.doc.pdf",width = 12,height = 136.5,limitsize = FALSE,
       path = "/Users/zacrobinson/Documents/For Publication PTF Meta/Plots")


