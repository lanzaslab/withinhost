title: "Sensitivity code for within host model of drug, E.coli, and P. multocida"

# run it to clear environment

rm(list=ls())

library("sensobol")
library("data.table")
library("ggplot2")
library(plyr)
library(tidyverse)
library(deSolve)
library("foreach")
library("parallel")
library("doParallel")
library("cowplot")
library("extrafont")

library(ggplot2)

budworm_fun <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    sigma = 0.28;
    Nmmax = 512861;
    d = 1.01;
    Nemax = 484577 ;
    lambda = 0.25;
    eta = 0.4;
    c=0.05; p=0.05; phi=0.14;
    signa <- data.frame(times = times, import = rep(0, length(times)))
    
    signa$import[signa$times==1]<-1000
    signa$import[signa$times==24]<-1000
    signa$import[signa$times==48]<-1000
    input <- approxfun(signa, rule = 2)
    
    dose= input(t)
    
    dS = -k*S+dose;
    dP = beta*k*S - (alpha+gamma)*P;
    dC = alpha*P - vartheta*C;
    dL = gamma*P - delta*L;
    dSe = sigma*(1-(Se+Re)/Nemax)*Se -d*((C/(C+Cs50))*Se);
    dRe = (1-c)*sigma*(1-(Se+Re)/Nemax)*Re - d*((C/(C+Cr50))*Re);
    dSm = lambda*(1-(Sm+Rm)/Nmmax)*Sm -eta*((L/(L+Ls50))*Sm)-phi*Sm;
    dRm = (1-p)*lambda*(1-(Sm+Rm)/Nmmax)*Rm - eta*((L/(L+Lr50))*Rm)-phi*Rm;
    list(c(dS, dP, dC, dL, dSe,dRe,dSm,dRm))
    
  })
}
times <- seq(0, 150, 1)

x01=c(S=0, P=0, C=0, L=0, Se=100000, Re=1000, Sm=100000, Rm=1000)

parms=c(sigma = 0.28,  Nmmax = 512861, d = 1.01,   Nemax = 484577,lambda = 0.25,eta = 0.4,
        c=0.05, p=0.05, phi=0.14, k=0.1 , beta=0.0008 , alpha=0.16  ,gamma=0.16 ,vartheta= 0.53, delta=0.24 , Cs50=0.45 , Cr50=0.47, Ls50= 0.016, Lr50=0.05 )

out.data = as.data.frame(lsoda(x01, times, budworm_fun , parms))

data<-out.data%>% select( time)
data1<- out.data$Rm
data$Re <-data1
ggplot(data = data, aes(x = time))+geom_line(aes(y = data1))

# Model for 5mg/kg three times dose
N <-2000
params <-c( "k", "beta", "alpha",  "gamma", "vartheta", "delta", "Cs50", "Cr50","Ls50","Lr50")
order <- "first"
R <- 10^3
#type <- "norm"
conf <- 0.95
times <- seq(0, 150, 1)
timeOutput <- seq(5, 150, 5)


budworm_fun <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    sigma = 0.28;
    Nmmax = 512861.91;
    d = 1.01;
    Nemax = 484577.76 ;
    lambda = 0.25;
    eta = 0.4;
    c=0.05; p=0.05; phi=0.14;
    signa <- data.frame(times = times, import = rep(0, length(times)))
    
    signa$import[signa$times==1]<-1000
    signa$import[signa$times==24]<-1000
    signa$import[signa$times==48]<-1000
    input <- approxfun(signa, rule = 2)
    
    dose= input(t)
    
    dS = -k*S+dose;
    dP = beta*k*S - (alpha+gamma)*P;
    dC = alpha*P - vartheta*C;
    dL = gamma*P - delta*L;
    dSe = sigma*(1-(Se+Re)/Nemax)*Se -d*((C/(C+Cs50))*Se);
    dRe = (1-c)*sigma*(1-(Se+Re)/Nemax)*Re - d*((C/(C+Cr50))*Re);
    dSm = lambda*(1-(Sm+Rm)/Nmmax)*Sm -eta*((L/(L+Ls50))*Sm)-phi*Sm;
    dRm = (1-p)*lambda*(1-(Sm+Rm)/Nmmax)*Rm - eta*((L/(L+Lr50))*Rm)-phi*Rm;
    list(c(dS, dP, dC, dL, dSe,dRe,dSm,dRm))
    
  })
}

#0.09, 0.9  0.07, 0.7 0.2, 0.7

#Parameters distribution


mat <- sobol_matrices(N = N, params = params, order = order,type="LHS" )
mat[, "k"]<- qunif(mat[, "k"], 0.04, 0.9)
mat[, "beta"] <- qunif(mat[, "beta"], 0.0001, 0.005)
mat[, "alpha"] <- qunif(mat[, "alpha"], 0.01, 0.5)
mat[, "gamma"] <- qunif(mat[, "gamma"], 0.01, 0.5)
mat[, "vartheta"] <- qunif(mat[, "vartheta"], 0.05, 0.9)
mat[, "delta"] <- qunif(mat[, "delta"], 0.07, 0.7)
mat[, "Cs50"] <- qunif(mat[, "Cs50"], 0.2, 0.45)
mat[, "Cr50"] <- qunif(mat[, "Cr50"], 0.3, 0.9)
mat[, "Ls50"] <- qunif(mat[, "Ls50"], 0.006, 0.06)
mat[, "Lr50"] <- qunif(mat[, "Lr50"], 0.05, 0.13)

#ODE solutions with input the uniform distribution of the parameters

n.cores <- makeCluster(floor(detectCores() * 0.75))
registerDoParallel(n.cores)
y <- foreach(i = 1:nrow(mat), .combine = "rbind",
             .packages = "sensobol") %dopar% {
               sobol_ode(d = mat[i, ], times = times, timeOutput = timeOutput,
                         state = c(S=0, P=0, C=0, L=0,Se=100000,Re=1000,Sm=100000,Rm=1000), func = budworm_fun)
             }
stopCluster(n.cores)

full.dt <- data.table(y)

indices.dt <- melt(full.dt, measure.vars = c("S", "P", "C", "L","Se","Re","Sm","Rm"))
#print(indices.dt)


#plotting of uncertainty

p.dt_sm<-indices.dt %>% filter(variable %in% c("Sm"))
p.dt_rm<-indices.dt %>% filter(variable %in% c("Rm"))
p.dt_se<-indices.dt %>% filter(variable %in% c("Se"))
p.dt_re<-indices.dt %>% filter(variable %in% c("Re"))

ggp_sm<-p.dt_sm%>%group_by(time) %>% 
  summarize(median = quantile(value, p = .5), q1 = quantile(value, p = .25), q2=quantile(value, p = .75),q3=quantile(value,p=.05), q4=quantile(value,p=.15), q5=quantile(value,p=.35),q6=quantile(value,p=.45),q7=quantile(value,p=.55),q8=quantile(value,p=.65),q9=quantile(value,p=.75),q10=quantile(value,p=.95))

ggp_rm<-p.dt_rm%>%group_by(time) %>% 
  summarize(median = quantile(value, p = .5), q1 = quantile(value, p = .25), q2=quantile(value, p = .75),q3=quantile(value,p=.05), q4=quantile(value,p=.15), q5=quantile(value,p=.35),q6=quantile(value,p=.45),q7=quantile(value,p=.55),q8=quantile(value,p=.65),q9=quantile(value,p=.75),q10=quantile(value,p=.95))

ggp_se<-p.dt_se%>%group_by(time) %>% 
  summarize(median = quantile(value, p = .5), q1 = quantile(value, p = .25), q2=quantile(value, p = .75),q3=quantile(value,p=.05), q4=quantile(value,p=.15), q5=quantile(value,p=.35),q6=quantile(value,p=.45),q7=quantile(value,p=.55),q8=quantile(value,p=.65),q9=quantile(value,p=.75),q10=quantile(value,p=.95))

ggp_re<-p.dt_re%>%group_by(time) %>% 
  summarize(median = quantile(value, p = .5), q1 = quantile(value, p = .25), q2=quantile(value, p = .75),q3=quantile(value,p=.05), q4=quantile(value,p=.15), q5=quantile(value,p=.35),q6=quantile(value,p=.45),q7=quantile(value,p=.55),q8=quantile(value,p=.65),q9=quantile(value,p=.75),q10=quantile(value,p=.95))

prob_plot_sm<- ggplot(ggp_sm, aes(x=time))  + geom_ribbon(aes(ymin = q1, ymax=q5,fill = "(25-35)"))+  geom_ribbon(aes(ymin = q5, ymax=q6,fill = "(35-45)"))+ geom_ribbon(aes(ymin = q6, ymax=median,fill = "(45-50)"))+ geom_ribbon(aes(ymin = median, ymax=q7,fill = "(50-55)"))+ geom_ribbon(aes(ymin = q7, ymax=q8,fill = "(55-65)"))+ geom_ribbon(aes(ymin = q8, ymax=q9,fill = "(65-75)"))+ geom_line(aes(y=median,colour="Median"),size=1) + scale_colour_manual(name='',values="darkred")+scale_fill_manual(values=c( "(25-35)"="#fee2f0", "(35-45)"="#FFB6C1", "(45-50)"="#FF1493", "(50-55)"="#FF1493", "(55-65)"="#FFB6C1", "(65-75)"="#fee2f0"))+ labs(fill="Quantile")+theme_minimal()+theme_bw() + theme( panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.title = element_text(face="bold",size = 14, family="Times New Roman"))+ theme(axis.text.x = element_text(face="bold", size=14, family="Times New Roman"),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   axis.text.y = element_text(face="bold", size=14, family="Times New Roman")) +labs(title="(b)",  x = "Time (hrs)", y= "Sm" )+scale_y_continuous(labels=scales::comma_format(scale=1))+theme(plot.title = element_text(face="bold",size = 14, family="Times New Roman"))

prob_plot_rm<- ggplot(ggp_rm, aes(x=time))  + geom_ribbon(aes(ymin = q1, ymax=q5,fill = "(25-35)"))+  geom_ribbon(aes(ymin = q5, ymax=q6,fill = "(35-45)"))+ geom_ribbon(aes(ymin = q6, ymax=median,fill = "(45-50)"))+ geom_ribbon(aes(ymin = median, ymax=q7,fill = "(50-55)"))+ geom_ribbon(aes(ymin = q7, ymax=q8,fill = "(55-65)"))+ geom_ribbon(aes(ymin = q8, ymax=q9,fill = "(65-75)"))+ geom_line(aes(y=median,colour="Median"),size=1) + scale_colour_manual(name='',values="darkred")+scale_fill_manual(values=c( "(25-35)"="#fee2f0", "(35-45)"="#FFB6C1", "(45-50)"="#FF1493", "(50-55)"="#FF1493", "(55-65)"="#FFB6C1", "(65-75)"="#fee2f0"))+ labs(fill="Quantile")+theme_minimal()+theme_bw() + theme( panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.title = element_text(face="bold",size = 14, family="Times New Roman"))+ theme(axis.text.x = element_text(face="bold", size=14, family="Times New Roman"),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   axis.text.y = element_text(face="bold", size=14, family="Times New Roman")) +labs( title="(d)", x = "Time (hrs)", y= "Rm" )+scale_y_continuous(labels=scales::comma_format(scale=1))+theme(plot.title = element_text(face="bold",size = 14, family="Times New Roman"))

prob_plot_se<- ggplot(ggp_se, aes(x=time))  + geom_ribbon(aes(ymin = q1, ymax=q5,fill = "(25-35)"))+  geom_ribbon(aes(ymin = q5, ymax=q6,fill = "(35-45)"))+ geom_ribbon(aes(ymin = q6, ymax=median,fill = "(45-50)"))+ geom_ribbon(aes(ymin = median, ymax=q7,fill = "(50-55)"))+ geom_ribbon(aes(ymin = q7, ymax=q8,fill = "(55-65)"))+ geom_ribbon(aes(ymin = q8, ymax=q9,fill = "(65-75)"))+ geom_line(aes(y=median,colour="Median"),size=1) + scale_colour_manual(name='',values="darkred")+scale_fill_manual(values=c( "(25-35)"="#fee2f0", "(35-45)"="#FFB6C1", "(45-50)"="#FF1493", "(50-55)"="#FF1493", "(55-65)"="#FFB6C1", "(65-75)"="#fee2f0"))+ labs(fill="Quantile")+theme_minimal()+theme_bw() + theme( panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.title = element_text(face="bold",size = 14, family="Times New Roman"))+ theme(axis.text.x = element_text(face="bold", size=14, family="Times New Roman"),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   axis.text.y = element_text(face="bold", size=14, family="Times New Roman")) +labs(title="(a)",  x = "Time (hrs)", y= "Se" )+scale_y_continuous(labels=scales::comma_format(scale=1))+theme(plot.title = element_text(face="bold",size = 14, family="Times New Roman"))

prob_plot_re<- ggplot(ggp_re, aes(x=time))  + geom_ribbon(aes(ymin = q1, ymax=q5,fill = "(25-35)"))+  geom_ribbon(aes(ymin = q5, ymax=q6,fill = "(35-45)"))+ geom_ribbon(aes(ymin = q6, ymax=median,fill = "(45-50)"))+ geom_ribbon(aes(ymin = median, ymax=q7,fill = "(50-55)"))+ geom_ribbon(aes(ymin = q7, ymax=q8,fill = "(55-65)"))+ geom_ribbon(aes(ymin = q8, ymax=q9,fill = "(65-75)"))+ geom_line(aes(y=median,colour="Median"),size=1) + scale_colour_manual(name='',values="darkred")+scale_fill_manual(values=c( "(25-35)"="#fee2f0", "(35-45)"="#FFB6C1", "(45-50)"="#FF1493", "(50-55)"="#FF1493", "(55-65)"="#FFB6C1", "(65-75)"="#fee2f0"))+ labs(fill="Quantile")+theme_minimal()+theme_bw() + theme( panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.title = element_text(face="bold",size = 14, family="Times New Roman"))+ theme(axis.text.x = element_text(face="bold", size=14, family="Times New Roman"),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   axis.text.y = element_text(face="bold", size=14, family="Times New Roman")) +labs( title= "(c)", x = "Time (hrs)", y=  "Re" )+scale_y_continuous(labels=scales::comma_format(scale=1))+theme(plot.title = element_text(face="bold",size = 14, family="Times New Roman"))

prob_plot<- plot_grid(prob_plot_se, prob_plot_sm, prob_plot_re, prob_plot_rm ,
                      
                      ncol = 2, nrow = 2)

ggsave("dis_1000.png", plot = prob_plot, dpi = 300, width = 10, height = 8)



# Sobol indices 

ncpus <- floor(detectCores() * 0.75)
indices <- indices.dt[, sobol_indices(Y = value, N = N, params = params,
                                      order = order, boot = TRUE, first = "jansen", R = R,
                                      parallel = "multicore", ncpus = ncpus)$results, .(variable, time)]

p_re<-indices %>% filter(variable %in% c("Re")) 
p_rm<-indices %>% filter(variable %in% c("Rm")) 

Rm_T_re<-p_re %>% filter(sensitivity %in% c("Ti")) #total-order indices
Rm_S_re<-p_re %>% filter(sensitivity %in% c("Si")) #first-order indices

Rm_T_rm<-p_rm %>% filter(sensitivity %in% c("Ti")) #total-order indices
Rm_S_rm<-p_rm %>% filter(sensitivity %in% c("Si")) #first-order indices


#plotting of indices
#total-order indices

gg_re<-ggplot(Rm_T_re, aes(x=time)) + geom_line(aes(y=original, 
                                                    color = parameters, group = parameters),size=1)+
  scale_colour_manual(values=c("lightblue","black","blue","pink","cyan", "purple", "green","chocolate","orange","brown1"))+labs(title =  expression(paste("Total-order sensitivity analysis for resistant" , italic(  " E.coli"), " (Re)")),  x = "Time (hrs)", y= "Sobol indices")+ theme_minimal()+guides(color=guide_legend(override.aes=list(size=1)))+  theme_bw()+ theme( panel.grid.major = element_blank(),
                                                                                                                                                                                                                                                                                                                                                                                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.title = element_text(face="bold", size = 12, family="Times New Roman"))+ theme(axis.text.x = element_text(face="bold", size=14, family="Times New Roman"),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          axis.text.y = element_text(face="bold", size=14, family="Times New Roman")) +theme(plot.title = element_text(face="bold", size = 11, family="Times New Roman"))
#gg

#ggsave("sen_tot_re_1000.png", plot = gg, dpi = 300, width = 8, height = 8)

#first-order indices

gg_s_re<-ggplot(Rm_S_re, aes(x=time)) + geom_line(aes(y=original, 
                                                      color = parameters, group = parameters),size=1)+scale_colour_manual(values=c("lightblue","black","blue","pink","cyan", "purple", "green","chocolate","orange","brown1"))+labs(title = expression(paste("First-order sensitivity analysis for resistant" , italic(  " E.coli"), "(Re)")),  x = "Time (hrs)", y= "Sobol indices")+ theme_minimal()+guides(color=guide_legend(override.aes=list(size=1)))+  theme_bw()+ theme( panel.grid.major = element_blank(),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.title = element_text(face="bold", size = 12, family="Times New Roman"))+ theme(axis.text.x = element_text(face="bold", size=14, family="Times New Roman"),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            axis.text.y = element_text(face="bold", size=14, family="Times New Roman"))+theme(plot.title = element_text(face="bold", size = 11, family="Times New Roman")) 



gg_rm<-ggplot(Rm_T_rm, aes(x=time)) + geom_line(aes(y=original, 
                                                    color = parameters, group = parameters),size=1)+
  scale_colour_manual(values=c("lightblue","black","blue","pink","cyan", "purple", "green","chocolate","orange","brown1"))+labs(title =  expression(paste("Total-order sensitivity analysis for resistant" , italic(  " P.multocida"), " (Rm)")),  x = "Time (hrs)", y= "Sobol indices")+ theme_minimal()+guides(color=guide_legend(override.aes=list(size=1)))+  theme_bw()+ theme( panel.grid.major = element_blank(),
                                                                                                                                                                                                                                                                                                                                                                                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.title = element_text(face="bold", size = 12, family="Times New Roman"))+ theme(axis.text.x = element_text(face="bold", size=14, family="Times New Roman"),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               axis.text.y = element_text(face="bold", size=14, family="Times New Roman")) +theme(plot.title = element_text(face="bold", size = 11, family="Times New Roman"))
#gg

#ggsave("sen_tot_re_1000.png", plot = gg, dpi = 300, width = 8, height = 8)

#first-order indices

gg_s_rm<-ggplot(Rm_S_rm, aes(x=time)) + geom_line(aes(y=original, 
                                                      color = parameters, group = parameters),size=1)+scale_colour_manual(values=c("lightblue","black","blue","pink","cyan", "purple", "green","chocolate","orange","brown1"))+labs(title = expression(paste("First-order sensitivity analysis for resistant" , italic(  " P.multocida"), "(Rm)")),  x = "Time (hrs)", y= "Sobol indices")+ theme_minimal()+guides(color=guide_legend(override.aes=list(size=1)))+  theme_bw()+ theme( panel.grid.major = element_blank(),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.title = element_text(face="bold", size = 12, family="Times New Roman"))+ theme(axis.text.x = element_text(face="bold", size=14, family="Times New Roman"),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 axis.text.y = element_text(face="bold", size=14, family="Times New Roman"))+theme(plot.title = element_text(face="bold", size = 11, family="Times New Roman")) 







#gg_s
gg_re_1000<- plot_grid(gg_s_re, gg_re, labels = c("(a)", "(b)"),   ncol = 2, nrow = 1)
gg_rm_1000<- plot_grid(gg_s_rm, gg_rm, labels = c("(a)", "(b)"),  ncol = 2, nrow = 1)

ggsave("sen_re_1000.png", plot = gg_re_1000, dpi = 300, width = 10, height = 4)
ggsave("sen_rm_1000.png", plot = gg_rm_1000, dpi = 300, width = 10, height = 4)


#integrated area of parameter


k_data_tot<-Rm_T %>% filter(parameters %in% c("k"))
k_data_fir<-Rm_S %>% filter(parameters %in% c("k"))

x_1_tot <- k_data_tot$time
y_1_tot <- abs(k_data_tot$original)

x_1_fir <- k_data_fir$time
y_1_fir <- abs(k_data_fir$original)

interp_function_k_tot <- approxfun(x_1_tot, y_1_tot)
interp_function_k_fir <- approxfun(x_1_fir, y_1_fir)

# Calculate the integral of the interpolated function over a specified range
result_k_tot <- integrate(interp_function_k_tot, lower = 5, upper = 150)

result_k_fir <- integrate(interp_function_k_fir, lower = 5, upper = 150)

# Print the result
k_fir <-cat("The integral of the data is re:", result_k_fir$value, "\n")
k_tot <-cat("The integral of the data is re:", result_k_tot$value, "\n")



# Print an estimate of the absolute error in the result
cat("Estimated absolute error:", result_k_tot$abs.error, "\n")
cat("Estimated absolute error:", result_k_fir$abs.error, "\n")


#sobol_indices(Y=full.dt[["Rm"]], N=N, params=params, order = order, boot = TRUE, R=R)



# model for 7.5 mg/kg single time dose 

N <-2000
params <-c( "k", "beta", "alpha",  "gamma", "vartheta", "delta", "Cs50", "Cr50","Ls50","Lr50")
order <- "first"
R <- 10^3
#type <- "norm"
conf <- 0.95
times <- seq(0, 150, 1)
timeOutput <- seq(5, 150, 5)


budworm_fun1 <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    sigma = 0.28;
    Nmmax = 512861.91;
    d = 1.01;
    Nemax = 484577.76 ;
    lambda = 0.25;
    eta = 0.4;
    c=0.05; p=0.05; phi=0.14;
    signa <- data.frame(times = times, import = rep(0, length(times)))
    
    signa$import[signa$times==1]<-1500
    #signa$import[signa$times==24]<-1000
    #signa$import[signa$times==48]<-1000
    input <- approxfun(signa, rule = 2)
    
    dose= input(t)
    
    dS = -k*S+dose;
    dP = beta*k*S - (alpha+gamma)*P;
    dC = alpha*P - vartheta*C;
    dL = gamma*P - delta*L;
    dSe = sigma*(1-(Se+Re)/Nemax)*Se -d*((C/(C+Cs50))*Se);
    dRe = (1-c)*sigma*(1-(Se+Re)/Nemax)*Re - d*((C/(C+Cr50))*Re);
    dSm = lambda*(1-(Sm+Rm)/Nmmax)*Sm -eta*((L/(L+Ls50))*Sm)-phi*Sm;
    dRm = (1-p)*lambda*(1-(Sm+Rm)/Nmmax)*Rm - eta*((L/(L+Lr50))*Rm)-phi*Rm;
    list(c(dS, dP, dC, dL, dSe,dRe,dSm,dRm))
    
  })
}
mat <- sobol_matrices(N = N, params = params, order = order,type="LHS" )
mat[, "k"]<- qunif(mat[, "k"], 0.04, 0.9)
mat[, "beta"] <- qunif(mat[, "beta"], 0.0001, 0.005)
mat[, "alpha"] <- qunif(mat[, "alpha"], 0.01, 0.5)
mat[, "gamma"] <- qunif(mat[, "gamma"], 0.01, 0.5)
mat[, "vartheta"] <- qunif(mat[, "vartheta"], 0.05, 0.9)
mat[, "delta"] <- qunif(mat[, "delta"], 0.07, 0.7)
mat[, "Cs50"] <- qunif(mat[, "Cs50"], 0.2, 0.45)
mat[, "Cr50"] <- qunif(mat[, "Cr50"], 0.3, 0.9)
mat[, "Ls50"] <- qunif(mat[, "Ls50"], 0.006, 0.06)
mat[, "Lr50"] <- qunif(mat[, "Lr50"], 0.05, 0.13)

n.cores <- makeCluster(floor(detectCores() * 0.75))
registerDoParallel(n.cores)
y1 <- foreach(i = 1:nrow(mat), .combine = "rbind",
              .packages = "sensobol") %dopar% {
                sobol_ode(d = mat[i, ], times = times, timeOutput = timeOutput,
                          state = c(S=0, P=0, C=0, L=0,Se=100000,Re=1000,Sm=100000,Rm=1000), func = budworm_fun1)
              }
stopCluster(n.cores)

full.dt1 <- data.table(y1)

indices.dt1 <- melt(full.dt1, measure.vars = c("S", "P", "C", "L","Se","Re","Sm","Rm"))
print(indices.dt1)

#plotting uncertainty


p.dt1_rm<-indices.dt1 %>% filter(variable %in% c("Rm")) 
p.dt1_re<-indices.dt1 %>% filter(variable %in% c("Re")) 
p.dt1_sm<-indices.dt1 %>% filter(variable %in% c("Sm")) 
p.dt1_se<-indices.dt1 %>% filter(variable %in% c("Se")) 


ggp1_rm<-p.dt1_rm%>%group_by(time) %>% 
  summarize(median = quantile(value, p = .5), q1 = quantile(value, p = .25), q2=quantile(value, p = .75),q3=quantile(value,p=.05), q4=quantile(value,p=.15), q5=quantile(value,p=.35),q6=quantile(value,p=.45),q7=quantile(value,p=.55),q8=quantile(value,p=.65),q9=quantile(value,p=.75),q10=quantile(value,p=.95))

ggp1_re<-p.dt1_re%>%group_by(time) %>% 
  summarize(median = quantile(value, p = .5), q1 = quantile(value, p = .25), q2=quantile(value, p = .75),q3=quantile(value,p=.05), q4=quantile(value,p=.15), q5=quantile(value,p=.35),q6=quantile(value,p=.45),q7=quantile(value,p=.55),q8=quantile(value,p=.65),q9=quantile(value,p=.75),q10=quantile(value,p=.95))

ggp1_sm<-p.dt1_sm%>%group_by(time) %>% 
  summarize(median = quantile(value, p = .5), q1 = quantile(value, p = .25), q2=quantile(value, p = .75),q3=quantile(value,p=.05), q4=quantile(value,p=.15), q5=quantile(value,p=.35),q6=quantile(value,p=.45),q7=quantile(value,p=.55),q8=quantile(value,p=.65),q9=quantile(value,p=.75),q10=quantile(value,p=.95))

ggp1_se<-p.dt1_se%>%group_by(time) %>% 
  summarize(median = quantile(value, p = .5), q1 = quantile(value, p = .25), q2=quantile(value, p = .75),q3=quantile(value,p=.05), q4=quantile(value,p=.15), q5=quantile(value,p=.35),q6=quantile(value,p=.45),q7=quantile(value,p=.55),q8=quantile(value,p=.65),q9=quantile(value,p=.75),q10=quantile(value,p=.95))



prob_plot1_rm<- ggplot(ggp1_rm, aes(x=time))  + geom_ribbon(aes(ymin = q1, ymax=q5,fill = "(25-35)"))+  geom_ribbon(aes(ymin = q5, ymax=q6,fill = "(35-45)"))+ geom_ribbon(aes(ymin = q6, ymax=median,fill = "(45-50)"))+ geom_ribbon(aes(ymin = median, ymax=q7,fill = "(50-55)"))+ geom_ribbon(aes(ymin = q7, ymax=q8,fill = "(55-65)"))+ geom_ribbon(aes(ymin = q8, ymax=q9,fill = "(65-75)"))+ geom_line(aes(y=median,colour="Median"),size=1) + scale_colour_manual(name='',values="darkred")+scale_fill_manual(values=c( "(25-35)"="#fee2f0", "(35-45)"="#FFB6C1", "(45-50)"="#FF1493", "(50-55)"="#FF1493", "(55-65)"="#FFB6C1", "(65-75)"="#fee2f0"))+ labs(fill="Quantile")+theme_minimal()+theme_bw() + theme( panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.title = element_text(face="bold",size = 14, family="Times New Roman")) +labs(  title="(d)",x = "Time (hrs)", y= "Rm" )+ theme(axis.text.x = element_text(face="bold", size=14, family="Times New Roman"),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     axis.text.y = element_text(face="bold", size=14, family="Times New Roman"))+scale_y_continuous(labels=scales::comma_format(scale=1))+theme(plot.title = element_text(face="bold",size = 14, family="Times New Roman"))


prob_plot1_re<- ggplot(ggp1_re, aes(x=time))  + geom_ribbon(aes(ymin = q1, ymax=q5,fill = "(25-35)"))+  geom_ribbon(aes(ymin = q5, ymax=q6,fill = "(35-45)"))+ geom_ribbon(aes(ymin = q6, ymax=median,fill = "(45-50)"))+ geom_ribbon(aes(ymin = median, ymax=q7,fill = "(50-55)"))+ geom_ribbon(aes(ymin = q7, ymax=q8,fill = "(55-65)"))+ geom_ribbon(aes(ymin = q8, ymax=q9,fill = "(65-75)"))+ geom_line(aes(y=median,colour="Median"),size=1) + scale_colour_manual(name='',values="darkred")+scale_fill_manual(values=c( "(25-35)"="#fee2f0", "(35-45)"="#FFB6C1", "(45-50)"="#FF1493", "(50-55)"="#FF1493", "(55-65)"="#FFB6C1", "(65-75)"="#fee2f0"))+ labs(fill="Quantile")+theme_minimal()+theme_bw() + theme( panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.title = element_text(face="bold",size = 14, family="Times New Roman")) +labs(  title="(c)",x = "Time (hrs)", y=  "Re" )+ theme(axis.text.x = element_text(face="bold", size=14, family="Times New Roman"),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      axis.text.y = element_text(face="bold", size=14, family="Times New Roman"))+scale_y_continuous(labels=scales::comma_format(scale=1)) +theme(plot.title = element_text(face="bold",size = 14, family="Times New Roman"))

prob_plot1_sm<- ggplot(ggp1_sm, aes(x=time))  + geom_ribbon(aes(ymin = q1, ymax=q5,fill = "(25-35)"))+  geom_ribbon(aes(ymin = q5, ymax=q6,fill = "(35-45)"))+ geom_ribbon(aes(ymin = q6, ymax=median,fill = "(45-50)"))+ geom_ribbon(aes(ymin = median, ymax=q7,fill = "(50-55)"))+ geom_ribbon(aes(ymin = q7, ymax=q8,fill = "(55-65)"))+ geom_ribbon(aes(ymin = q8, ymax=q9,fill = "(65-75)"))+ geom_line(aes(y=median,colour="Median"),size=1) + scale_colour_manual(name='',values="darkred")+scale_fill_manual(values=c( "(25-35)"="#fee2f0", "(35-45)"="#FFB6C1", "(45-50)"="#FF1493", "(50-55)"="#FF1493", "(55-65)"="#FFB6C1", "(65-75)"="#fee2f0"))+ labs(fill="Quantile")+theme_minimal()+theme_bw() + theme( panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.title = element_text(face="bold",size = 14)) +labs(  title="(b)",x = "Time (hrs)", y=  "Sm" )+ theme(axis.text.x = element_text(face="bold", size=14, family="Times New Roman"),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            axis.text.y = element_text(face="bold", size=14, family="Times New Roman"))+scale_y_continuous(labels=scales::comma_format(scale=1))+theme(plot.title = element_text(face="bold",size = 14, family="Times New Roman"))

prob_plot1_se<- ggplot(ggp1_se, aes(x=time))  + geom_ribbon(aes(ymin = q1, ymax=q5,fill = "(25-35)"))+  geom_ribbon(aes(ymin = q5, ymax=q6,fill = "(35-45)"))+ geom_ribbon(aes(ymin = q6, ymax=median,fill = "(45-50)"))+ geom_ribbon(aes(ymin = median, ymax=q7,fill = "(50-55)"))+ geom_ribbon(aes(ymin = q7, ymax=q8,fill = "(55-65)"))+ geom_ribbon(aes(ymin = q8, ymax=q9,fill = "(65-75)"))+ geom_line(aes(y=median,colour="Median"),size=1) + scale_colour_manual(name='',values="darkred")+scale_fill_manual(values=c( "(25-35)"="#fee2f0", "(35-45)"="#FFB6C1", "(45-50)"="#FF1493", "(50-55)"="#FF1493", "(55-65)"="#FFB6C1", "(65-75)"="#fee2f0"))+ labs(fill="Quantile")+theme_minimal()+theme_bw() + theme( panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.title = element_text(face="bold",size = 14)) +labs( title="(a)", x = "Time (hrs)", y= "Se" )+ theme(axis.text.x = element_text(face="bold", size=14, family="Times New Roman"),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           axis.text.y = element_text(face="bold", size=14, family="Times New Roman"))+scale_y_continuous(labels=scales::comma_format(scale=1))+theme(plot.title = element_text(face="bold",size = 14, family="Times New Roman"))

prob_plot1<- plot_grid(prob_plot1_se, prob_plot1_sm, prob_plot1_re, prob_plot1_rm ,
                       
                       ncol = 2, nrow = 2)

ggsave("dis_1500.png", plot = prob_plot1, dpi = 300, width = 10, height = 8)



ncpus <- floor(detectCores() * 0.75)
indices1 <- indices.dt1[, sobol_indices(Y = value, N = N, params = params,
                                        order = order, boot = TRUE, first = "jansen", R = R,
                                        parallel = "multicore", ncpus = ncpus)$results, .(variable, time)]


p1_re<-indices1 %>% filter(variable %in% c("Re")) 
p1_rm<-indices1 %>% filter(variable %in% c("Rm"))



Rm_T1_re<-p1_re %>% filter(sensitivity %in% c("Ti"))
Rm_S1_re<-p1_re %>% filter(sensitivity %in% c("Si"))

Rm_T1_rm<-p1_rm %>% filter(sensitivity %in% c("Ti"))
Rm_S1_rm<-p1_rm %>% filter(sensitivity %in% c("Si")) 


#plotting indices

gg1_re<-ggplot(Rm_T1_re, aes(x=time)) + geom_line(aes(y=original, 
                                                      color = parameters, group = parameters),size=1)+scale_colour_manual(values=c("lightblue","black","blue","pink","cyan", "purple", "green","chocolate","orange","brown1"))+labs(title =  expression(paste("Total-order sensitivity analysis for resistant" , italic(  " E.coli"), "(Re)")),  x = "Time (hrs)", y= "Sobol indices")+ theme_minimal()+guides(color=guide_legend(override.aes=list(size=1.5)))+  theme_bw()+ theme( panel.grid.major = element_blank(),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ theme(axis.text.x = element_text(face="bold", size=12, family="Times New Roman"),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            axis.text.y = element_text(face="bold", size=12))+theme(axis.title = element_text(face="bold", size = 14, family="Times New Roman")) +theme(plot.title = element_text(face="bold", size = 11, family="Times New Roman"))


gg_s1_re<-ggplot(Rm_S1_re, aes(x=time)) + geom_line(aes(y=original, 
                                                        color = parameters, group = parameters),size=1)+scale_colour_manual(values=c("lightblue","black","blue","pink","cyan", "purple", "green","chocolate","orange","brown1"))+labs(title = expression(paste("First-order sensitivity analysis for resistant" , italic(  " E.coli"), "(Re)")),  x = "Time (hrs)", y= "Sobol indices")+ theme_minimal()+guides(color=guide_legend(override.aes=list(size=1.5)))+  theme_bw()+ theme( panel.grid.major = element_blank(),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.title = element_text(face="bold", size = 14, family="Times New Roman"))+ theme(axis.text.x = element_text(face="bold", size=12, family="Times New Roman"),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                axis.text.y = element_text(face="bold", size=12, family="Times New Roman"))+theme(plot.title = element_text(face="bold", size = 11, family="Times New Roman")) 


gg1_rm<-ggplot(Rm_T1_rm, aes(x=time)) + geom_line(aes(y=original, 
                                                      color = parameters, group = parameters),size=1)+scale_colour_manual(values=c("lightblue","black","blue","pink","cyan", "purple", "green","chocolate","orange","brown1"))+labs(title =  expression(paste("Total-order sensitivity analysis for resistant" , italic(  " P.multocida"), "(Rm)")),  x = "Time (hrs)", y= "Sobol indices")+ theme_minimal()+guides(color=guide_legend(override.aes=list(size=1.5)))+  theme_bw()+ theme( panel.grid.major = element_blank(),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ theme(axis.text.x = element_text(face="bold", size=12, family="Times New Roman"),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 axis.text.y = element_text(face="bold", size=12, family="Times New Roman"))+theme(axis.title = element_text(face="bold", size = 14, family="Times New Roman")) +theme(plot.title = element_text(face="bold", size = 11, family="Times New Roman"))


gg_s1_rm<-ggplot(Rm_S1_rm, aes(x=time)) + geom_line(aes(y=original, 
                                                        color = parameters, group = parameters),size=1)+scale_colour_manual(values=c("lightblue","black","blue","pink","cyan", "purple", "green","chocolate","orange","brown1"))+labs(title = expression(paste("First-order sensitivity analysis for resistant" , italic(  " P.multocida"), "(Rm)")),  x = "Time (hrs)", y= "Sobol indices")+ theme_minimal()+guides(color=guide_legend(override.aes=list(size=1.5)))+  theme_bw()+ theme( panel.grid.major = element_blank(),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.title = element_text(face="bold", size = 14, family="Times New Roman"))+ theme(axis.text.x = element_text(face="bold", size=12, family="Times New Roman"),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     axis.text.y = element_text(face="bold", size=12, family="Times New Roman"))+theme(plot.title = element_text(face="bold", size = 11, family="Times New Roman"))

#gg_s
gg_re_1500<- plot_grid(gg_s1_re, gg1_re, labels = c("(a)", "(b)"),  ncol = 2, nrow = 1)
gg_rm_1500<- plot_grid(gg_s1_rm, gg1_rm, labels = c("(a)", "(b)"),  ncol = 2, nrow = 1)

ggsave("sen_re_1500.png", plot = gg_re_1500, dpi = 300, width = 10, height = 4)
ggsave("sen_rm_1500.png", plot = gg_rm_1500, dpi = 300, width = 10, height = 4)


#integrated area of parameter
k1_data_tot<-Rm_T1 %>% filter(parameters %in% c("k"))
k1_data_fir<-Rm_S1 %>% filter(parameters %in% c("k"))

x_1_tot <- k1_data_tot$time
y_1_tot <-abs( k1_data_tot$original)

x_1_fir <- k1_data_fir$time
y_1_fir <-abs( k1_data_fir$original)

interp_function_k_tot <- approxfun(x_1_tot, y_1_tot)
interp_function_k_fir <- approxfun(x_1_fir, y_1_fir)

# Calculate the integral of the interpolated function over a specified range
result_k_tot <- integrate(interp_function_k_tot, lower = 5, upper = 150)

result_k_fir <- integrate(interp_function_k_fir, lower = 5, upper = 150)

# Print the result
k_fir <-cat("The integral of the data is_fir:", result_k_fir$value, "\n")
k_tot <-cat("The integral of the data is_tot:", result_k_tot$value, "\n")



# Print an estimate of the absolute error in the result
cat("Estimated absolute error:", result_k_tot$abs.error, "\n")
cat("Estimated absolute error:", result_k_fir$abs.error, "\n")




# plotting when drug is 12.5 mg/kg

library("sensobol")
library("data.table")
library("ggplot2")
N <-2000
params <-c( "k", "beta", "alpha",  "gamma", "vartheta", "delta", "Cs50", "Cr50","Ls50","Lr50")
order <- "first"
R <- 10^3
#type <- "norm"
conf <- 0.95
times <- seq(0, 150, 1)
timeOutput <- seq(5, 150, 5)


budworm_fun2 <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    sigma = 0.28;
    Nmmax = 512861.91;
    d = 1.01;
    Nemax = 484577.76 ;
    lambda = 0.25;
    eta = 0.4;
    c=0.05; p=0.05; phi=0.14;
    signa <- data.frame(times = times, import = rep(0, length(times)))
    
    signa$import[signa$times==1]<-2500
    #signa$import[signa$times==24]<-1000
    #signa$import[signa$times==48]<-1000
    input <- approxfun(signa, rule = 2)
    
    dose= input(t)
    
    dS = -k*S+dose;
    dP = beta*k*S - (alpha+gamma)*P;
    dC = alpha*P - vartheta*C;
    dL = gamma*P - delta*L;
    dSe = sigma*(1-(Se+Re)/Nemax)*Se -d*((C/(C+Cs50))*Se);
    dRe = (1-c)*sigma*(1-(Se+Re)/Nemax)*Re - d*((C/(C+Cr50))*Re);
    dSm = lambda*(1-(Sm+Rm)/Nmmax)*Sm -eta*((L/(L+Ls50))*Sm)-phi*Sm;
    dRm = (1-p)*lambda*(1-(Sm+Rm)/Nmmax)*Rm - eta*((L/(L+Lr50))*Rm)-phi*Rm;
    list(c(dS, dP, dC, dL, dSe,dRe,dSm,dRm))
    
  })
}
mat <- sobol_matrices(N = N, params = params, order = order,type="LHS" )
mat[, "k"]<- qunif(mat[, "k"], 0.04, 0.9)
mat[, "beta"] <- qunif(mat[, "beta"], 0.0001, 0.005)
mat[, "alpha"] <- qunif(mat[, "alpha"], 0.01, 0.5)
mat[, "gamma"] <- qunif(mat[, "gamma"], 0.01, 0.5)
mat[, "vartheta"] <- qunif(mat[, "vartheta"], 0.05, 0.9)
mat[, "delta"] <- qunif(mat[, "delta"], 0.07, 0.7)
mat[, "Cs50"] <- qunif(mat[, "Cs50"], 0.2, 0.45)
mat[, "Cr50"] <- qunif(mat[, "Cr50"], 0.3, 0.9)
mat[, "Ls50"] <- qunif(mat[, "Ls50"], 0.006, 0.06)
mat[, "Lr50"] <- qunif(mat[, "Lr50"], 0.05, 0.13)


library("foreach")
library("parallel")
library("doParallel")
library(plyr)
library(tidyverse)
library(deSolve)

n.cores <- makeCluster(floor(detectCores() * 0.75))
registerDoParallel(n.cores)
y2 <- foreach(i = 1:nrow(mat), .combine = "rbind",
              .packages = "sensobol") %dopar% {
                sobol_ode(d = mat[i, ], times = times, timeOutput = timeOutput,
                          state = c(S=0, P=0, C=0, L=0,Se=100000,Re=1000,Sm=100000,Rm=1000), func = budworm_fun2)
              }
stopCluster(n.cores)

full.dt2 <- data.table(y2)
#full.dt2 <- data.table(full.dt2[full.dt2$time==125,])

indices.dt2 <- melt(full.dt2, measure.vars = c("S", "P", "C", "L","Se","Re","Sm","Rm"))
print(indices.dt2)


#plotting uncertainty
p.dt2_sm<-indices.dt2 %>% filter(variable %in% c("Sm")) 
p.dt2_rm<-indices.dt2 %>% filter(variable %in% c("Rm"))
p.dt2_se<-indices.dt2 %>% filter(variable %in% c("Se")) 
p.dt2_re<-indices.dt2 %>% filter(variable %in% c("Re")) 

#max(p.dt2$value)
#min(p.dt2$value)

ggp2_sm<-p.dt2_sm%>%group_by(time) %>% 
  summarize(median = quantile(value, p = .5), q1 = quantile(value, p = .25), q2=quantile(value, p = .75),q3=quantile(value,p=.05), q4=quantile(value,p=.15), q5=quantile(value,p=.35),q6=quantile(value,p=.45),q7=quantile(value,p=.55),q8=quantile(value,p=.65),q9=quantile(value,p=.75),q10=quantile(value,p=.95))

ggp2_rm<-p.dt2_rm%>%group_by(time) %>% 
  summarize(median = quantile(value, p = .5), q1 = quantile(value, p = .25), q2=quantile(value, p = .75),q3=quantile(value,p=.05), q4=quantile(value,p=.15), q5=quantile(value,p=.35),q6=quantile(value,p=.45),q7=quantile(value,p=.55),q8=quantile(value,p=.65),q9=quantile(value,p=.75),q10=quantile(value,p=.95))

ggp2_se<-p.dt2_se%>%group_by(time) %>% 
  summarize(median = quantile(value, p = .5), q1 = quantile(value, p = .25), q2=quantile(value, p = .75),q3=quantile(value,p=.05), q4=quantile(value,p=.15), q5=quantile(value,p=.35),q6=quantile(value,p=.45),q7=quantile(value,p=.55),q8=quantile(value,p=.65),q9=quantile(value,p=.75),q10=quantile(value,p=.95))

ggp2_re<-p.dt2_re%>%group_by(time) %>% 
  summarize(median = quantile(value, p = .5), q1 = quantile(value, p = .25), q2=quantile(value, p = .75),q3=quantile(value,p=.05), q4=quantile(value,p=.15), q5=quantile(value,p=.35),q6=quantile(value,p=.45),q7=quantile(value,p=.55),q8=quantile(value,p=.65),q9=quantile(value,p=.75),q10=quantile(value,p=.95))


prob_plot2_sm<- ggplot(ggp2_sm, aes(x=time))  + geom_ribbon(aes(ymin = q1, ymax=q5,fill = "(25-35)"))+  geom_ribbon(aes(ymin = q5, ymax=q6,fill = "(35-45)"))+ geom_ribbon(aes(ymin = q6, ymax=median,fill = "(45-50)"))+ geom_ribbon(aes(ymin = median, ymax=q7,fill = "(50-55)"))+ geom_ribbon(aes(ymin = q7, ymax=q8,fill = "(55-65)"))+ geom_ribbon(aes(ymin = q8, ymax=q9,fill = "(65-75)"))+ geom_line(aes(y=median, colour="Median"),size=1) + scale_colour_manual(name='',values="darkred")+scale_fill_manual(values=c( "(25-35)"="#fee2f0", "(35-45)"="#FFB6C1", "(45-50)"="#FF1493", "(50-55)"="#FF1493", "(55-65)"="#FFB6C1", "(65-75)"="#fee2f0"))+ labs(fill="Quantile")+theme_minimal()+theme_bw() + theme( panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.title = element_text(face="bold", size = 14, family="Times New Roman"))+labs(title="(b)", x = "Time (hrs)", y= " Sm")+ theme(axis.text.x = element_text(face="bold", size=14, family="Times New Roman"),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     axis.text.y = element_text(face="bold", size=14, family="Times New Roman")) +scale_y_continuous(labels=scales::comma_format(scale=1))+theme(plot.title = element_text(face="bold",size = 14, family="Times New Roman"))


prob_plot2_rm<- ggplot(ggp2_rm, aes(x=time))  + geom_ribbon(aes(ymin = q1, ymax=q5,fill = "(25-35)"))+  geom_ribbon(aes(ymin = q5, ymax=q6,fill = "(35-45)"))+ geom_ribbon(aes(ymin = q6, ymax=median,fill = "(45-50)"))+ geom_ribbon(aes(ymin = median, ymax=q7,fill = "(50-55)"))+ geom_ribbon(aes(ymin = q7, ymax=q8,fill = "(55-65)"))+ geom_ribbon(aes(ymin = q8, ymax=q9,fill = "(65-75)"))+ geom_line(aes(y=median, colour="Median"),size=1) + scale_colour_manual(name='',values="darkred")+scale_fill_manual(values=c( "(25-35)"="#fee2f0", "(35-45)"="#FFB6C1", "(45-50)"="#FF1493", "(50-55)"="#FF1493", "(55-65)"="#FFB6C1", "(65-75)"="#fee2f0"))+ labs(fill="Quantile")+theme_minimal()+theme_bw() + theme( panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.title = element_text(face="bold", size = 14, family="Times New Roman"))+labs(title="(d)", x = "Time (hrs)", y= " Rm")+ theme(axis.text.x = element_text(face="bold", size=14, family="Times New Roman"),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     axis.text.y = element_text(face="bold", size=14, family="Times New Roman")) +scale_y_continuous(labels=scales::comma_format(scale=1))+theme(plot.title = element_text(face="bold",size = 14, family="Times New Roman"))




prob_plot2_se<- ggplot(ggp2_se, aes(x=time))  + geom_ribbon(aes(ymin = q1, ymax=q5,fill = "(25-35)"))+  geom_ribbon(aes(ymin = q5, ymax=q6,fill = "(35-45)"))+ geom_ribbon(aes(ymin = q6, ymax=median,fill = "(45-50)"))+ geom_ribbon(aes(ymin = median, ymax=q7,fill = "(50-55)"))+ geom_ribbon(aes(ymin = q7, ymax=q8,fill = "(55-65)"))+ geom_ribbon(aes(ymin = q8, ymax=q9,fill = "(65-75)"))+ geom_line(aes(y=median, colour="Median"),size=1) + scale_colour_manual(name='',values="darkred")+scale_fill_manual(values=c( "(25-35)"="#fee2f0", "(35-45)"="#FFB6C1", "(45-50)"="#FF1493", "(50-55)"="#FF1493", "(55-65)"="#FFB6C1", "(65-75)"="#fee2f0"))+ labs(fill="Quantile")+theme_minimal()+theme_bw() + theme( panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.title = element_text(face="bold", size = 14, family="Times New Roman"))+labs(title="(a)", x = "Time (hrs)", y= " Se")+ theme(axis.text.x = element_text(face="bold", size=14, family="Times New Roman"),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     axis.text.y = element_text(face="bold", size=14, family="Times New Roman")) +scale_y_continuous(labels=scales::comma_format(scale=1))+theme(plot.title = element_text(face="bold",size = 14, family="Times New Roman"))


prob_plot2_re<- ggplot(ggp2_re, aes(x=time))  + geom_ribbon(aes(ymin = q1, ymax=q5,fill = "(25-35)"))+  geom_ribbon(aes(ymin = q5, ymax=q6,fill = "(35-45)"))+ geom_ribbon(aes(ymin = q6, ymax=median,fill = "(45-50)"))+ geom_ribbon(aes(ymin = median, ymax=q7,fill = "(50-55)"))+ geom_ribbon(aes(ymin = q7, ymax=q8,fill = "(55-65)"))+ geom_ribbon(aes(ymin = q8, ymax=q9,fill = "(65-75)"))+ geom_line(aes(y=median, colour="Median"),size=1) + scale_colour_manual(name='',values="darkred")+scale_fill_manual(values=c( "(25-35)"="#fee2f0", "(35-45)"="#FFB6C1", "(45-50)"="#FF1493", "(50-55)"="#FF1493", "(55-65)"="#FFB6C1", "(65-75)"="#fee2f0"))+ labs(fill="Quantile")+theme_minimal()+theme_bw() + theme( panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.title = element_text(face="bold", size = 14, family="Times New Roman"))+labs(title="(c)", x = "Time (hrs)", y= " Re")+ theme(axis.text.x = element_text(face="bold", size=14, family="Times New Roman"),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     axis.text.y = element_text(face="bold", size=14, family="Times New Roman")) +scale_y_continuous(labels=scales::comma_format(scale=1))+theme(plot.title = element_text(face="bold",size = 14, family="Times New Roman"))




prob_plot2<- plot_grid(prob_plot2_se, prob_plot2_sm, prob_plot2_re, prob_plot2_rm ,
                       ncol = 2, nrow = 2)

ggsave("dis_2500.png", plot = prob_plot2, dpi = 300, width = 10, height = 8)



ncpus <- floor(detectCores() * 0.75)
indices2 <- indices.dt2[, sobol_indices(Y = value, N = N, params = params,
                                        order = order, boot = TRUE, first = "jansen", R = R,
                                        parallel = "multicore", ncpus = ncpus)$results, .(variable, time)]


 
p2_re<-indices2 %>% filter(variable %in% c("Re")) 
p2_rm<-indices2 %>% filter(variable %in% c("Rm")) 

Rm_T2_re<-p2_re %>% filter(sensitivity %in% c("Ti"))
Rm_S2_re<-p2_re %>% filter(sensitivity %in% c("Si"))

Rm_T2_rm<-p2_rm %>% filter(sensitivity %in% c("Ti"))
Rm_S2_rm<-p2_rm %>% filter(sensitivity %in% c("Si"))


#plotting indices 
gg2_re<-ggplot(Rm_T2_re, aes(x=time)) + geom_line(aes(y=original, 
                                                      color = parameters, group = parameters),size=1)+scale_colour_manual(values=c("lightblue","black","blue","pink","cyan", "purple", "green","chocolate","orange","brown1"))+labs(title =  expression(paste("Total-order sensitivity analysis for resistant" , italic(  " E.coli"), " (Re)")),  x = "Time (hrs)", y= "Sobol indices")+ theme_minimal()+guides(color=guide_legend(override.aes=list(size=1.5)))+  theme_bw()+ theme( panel.grid.major = element_blank(),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ theme(axis.text.x = element_text(face="bold", size=12, family="Times New Roman"),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             axis.text.y = element_text(face="bold", size=12, family="Times New Roman")) + theme(axis.title = element_text(face="bold", size = 14, family="Times New Roman"))+theme(plot.title = element_text(face="bold", size = 11, family="Times New Roman"))


gg2_s_re<-ggplot(Rm_S2_re, aes(x=time)) + geom_line(aes(y=original, 
                                                        color = parameters, group = parameters),size=1)+scale_colour_manual(values=c("lightblue","black","blue","pink","cyan", "purple", "green","chocolate","orange","brown1"))+labs(title =  expression(paste("First-order sensitivity analysis for resistant" , italic(  " E.coli"), " (Re)")),  x = "Time (hrs)", y= "Sobol indices")+ theme_minimal()+guides(color=guide_legend(override.aes=list(size=1.5)))+  theme_bw()+ theme( panel.grid.major = element_blank(),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ theme(axis.text.x = element_text(face="bold", size=12, family="Times New Roman"),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               axis.text.y = element_text(face="bold", size=12, family="Times New Roman")) + theme(axis.title = element_text(face="bold", size = 14, family="Times New Roman"))+theme(plot.title = element_text(face="bold", size = 11, family="Times New Roman"))


gg2_rm<-ggplot(Rm_T2_rm, aes(x=time)) + geom_line(aes(y=original, 
                                                      color = parameters, group = parameters),size=1)+scale_colour_manual(values=c("lightblue","black","blue","pink","cyan", "purple", "green","chocolate","orange","brown1"))+labs(title =  expression(paste("Total-order sensitivity analysis for resistant" , italic(  " P.multocida"), " (Rm)")),  x = "Time (hrs)", y= "Sobol indices")+ theme_minimal()+guides(color=guide_legend(override.aes=list(size=1.5)))+  theme_bw()+ theme( panel.grid.major = element_blank(),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ theme(axis.text.x = element_text(face="bold", size=12, family="Times New Roman"),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  axis.text.y = element_text(face="bold", size=12, family="Times New Roman")) + theme(axis.title = element_text(face="bold", size = 14, family="Times New Roman"))+theme(plot.title = element_text(face="bold", size = 11, family="Times New Roman"))


gg2_s_rm<-ggplot(Rm_S2_rm, aes(x=time)) + geom_line(aes(y=original, 
                                                        color = parameters, group = parameters),size=1)+scale_colour_manual(values=c("lightblue","black","blue","pink","cyan", "purple", "green","chocolate","orange","brown1"))+labs(title =  expression(paste("First-order sensitivity analysis for resistant" , italic(  " P.multocida"), " (Rm)")),  x = "Time (hrs)", y= "Sobol indices")+ theme_minimal()+guides(color=guide_legend(override.aes=list(size=1.5)))+  theme_bw()+ theme( panel.grid.major = element_blank(),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ theme(axis.text.x = element_text(face="bold", size=12, family="Times New Roman"),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    axis.text.y = element_text(face="bold", size=12, family="Times New Roman")) + theme(axis.title = element_text(face="bold", size = 14, family="Times New Roman"))+theme(plot.title = element_text(face="bold", size = 11, family="Times New Roman"))

gg_re_2500<- plot_grid(gg2_s_re, gg2_re, labels = c("(a)", "(b)"),  ncol = 2, nrow = 1)
gg_rm_2500<- plot_grid(gg2_s_rm, gg2_rm, labels = c("(a)", "(b)"),  ncol = 2, nrow = 1)

ggsave("sen_re_2500.png", plot = gg_re_2500, dpi = 300, width = 10, height = 4)
ggsave("sen_rm_2500.png", plot = gg_rm_2500, dpi = 300, width = 10, height = 4) 


#integrated area of parameter
k_data_tot<-Rm_T2 %>% filter(parameters %in% c("k"))
k_data_fir<-Rm_S2 %>% filter(parameters %in% c("k"))

x_1_tot <- k_data_tot$time
y_1_tot <- abs( k_data_tot$original)

x_1_fir <- k_data_fir$time
y_1_fir <- abs(k_data_fir$original)
data_list<-as.data.frame(y_1_fir)

interp_function_k_tot <- approxfun(x_1_tot, y_1_tot)
interp_function_k_fir <- approxfun(x_1_fir, y_1_fir)

# Calculate the integral of the interpolated function over a specified range
result_k_tot <- integrate(interp_function_k_tot, lower = 5, upper = 150)

result_k_fir <- integrate(interp_function_k_fir, lower = 5, upper = 150)

# Print the resultlist_data<-as.data.frame(result_k_fir)
k_fir <-cat("The integral of the data is_fir:", result_k_fir$value, "\n")
k_tot <-cat("The integral of the data is_tot:", result_k_tot$value, "\n")



# Print an estimate of the absolute error in the result
cat("Estimated absolute error:", result_k_tot$abs.error, "\n")
cat("Estimated absolute error:", result_k_fir$abs.error, "\n")

