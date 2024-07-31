# Single high dose (12.5mg/kg bwt)
rm(list=ls())
library(tidyverse)
library(deSolve)
library(data.table)

library(ggplot2)
install.packages("extrafont")
library(extrafont)

install.packages("patchwork")
library(patchwork)


library(gridExtra)

# 12.5 mg/kg-once

times <- seq(0, 150, by = 1)
signa <- data.frame(times = times, import = rep(0, length(times)))

signa$import[signa$times==1]<-2500
input <- approxfun(signa, rule = 2)
input(c(0,1,1.95833333,2,2.04166667,3,24,125))

optimalmodel = function(t,x,parms){
  with(as.list(c(parms,x)), {
    
    dose= input(t) 
    
    dS = -k*S+dose;
    dP = beta*k*S - (alpha+gamma)*P;
    dC = alpha*P - vartheta*C;
    dL = gamma*P - delta*L;
    dSe = sigma*(1-(Se+Re)/Nemax)*Se -d*((C/(C+Cs50))*Se);
    dRe = (1-c)*sigma*(1-(Se+Re)/Nemax)*Re - d*((C/(C+Cr50))*Re);
    dSm = lambda*(1-(Sm+Rm)/Nmmax)*Sm -eta*((L/(L+Ls50))*Sm)-phi*Sm;
    dRm = (1-p)*lambda*(1-(Sm+Rm)/Nmmax)*Rm - eta*((L/(L+Lr50))*Rm)-phi*Rm;
    dY=c(dS, dP, dC, dL, dSe,dRe,dSm,dRm);
    list(dY, signal= dose)
  }
  )
}

parms=c(k=0.44, beta= 0.0016, alpha= 0.16, vartheta= 0.53, sigma=0.28, Nemax=484577,
        Cs50= 0.25, Cr50= 0.471, d= 1.01, gamma=0.16 , lambda= 0.25, Nmmax=512861 ,
        Ls50=0.016 , Lr50= 0.125,delta= 0.24, eta=0.40, c=0.05, p=0.05, phi=0.14)
xstart= c(S=0,P=0,C=0,L=0,Se=400000,Re=1000,Sm=40000,Rm=1000)

Single_high = as.data.frame ( lsoda(xstart, times, optimalmodel, parms))



# 5 mg/kg-three

times <- seq(0, 150, by = 1)
signa <- data.frame(times = times, import = rep(0, length(times)))

signa$import[signa$times==1]<-1000
signa$import[signa$times==25]<-1000
signa$import[signa$times==49]<-1000
input <- approxfun(signa, rule = 2)
input(c(0,1,1.95833333,2,2.04166667,3,24,125))

optimalmodel = function(t,x,parms){
  with(as.list(c(parms,x)), {
    
    dose= input(t) 
    
    dS = -k*S+dose;
    dP = beta*k*S - (alpha+gamma)*P;
    dC = alpha*P - vartheta*C;
    dL = gamma*P - delta*L;
    dSe = sigma*(1-(Se+Re)/Nemax)*Se -d*((C/(C+Cs50))*Se);
    dRe = (1-c)*sigma*(1-(Se+Re)/Nemax)*Re - d*((C/(C+Cr50))*Re);
    dSm = lambda*(1-(Sm+Rm)/Nmmax)*Sm -eta*((L/(L+Ls50))*Sm)-phi*Sm;
    dRm = (1-p)*lambda*(1-(Sm+Rm)/Nmmax)*Rm - eta*((L/(L+Lr50))*Rm)-phi*Rm;
    dY=c(dS, dP, dC, dL, dSe,dRe,dSm,dRm);
    list(dY, signal= dose)
  }
  )
}

parms=c(k=0.44, beta= 0.0016, alpha= 0.16, vartheta= 0.53, sigma=0.28, Nemax=484577,
        Cs50= 0.25, Cr50= 0.471, d= 1.01, gamma=0.16 , lambda= 0.25, Nmmax=512861 ,
        Ls50=0.016 , Lr50= 0.125,delta= 0.24, eta=0.40, c=0.05, p=0.05, phi=0.14)
xstart= c(S=0,P=0,C=0,L=0,Se=400000,Re=1000,Sm=40000,Rm=1000)

Multiple_low = as.data.frame ( lsoda(xstart, times, optimalmodel, parms))




# 7.5 mg/kg-once

times <- seq(0, 150, by = 1)
signa <- data.frame(times = times, import = rep(0, length(times)))

signa$import[signa$times==1]<-1500
input <- approxfun(signa, rule = 2)
input(c(0,1,1.95833333,2,2.04166667,3,24,125))

optimalmodel = function(t,x,parms){
  with(as.list(c(parms,x)), {
    
    dose= input(t) 
    
    dS = -k*S+dose;
    dP = beta*k*S - (alpha+gamma)*P;
    dC = alpha*P - vartheta*C;
    dL = gamma*P - delta*L;
    dSe = sigma*(1-(Se+Re)/Nemax)*Se -d*((C/(C+Cs50))*Se);
    dRe = (1-c)*sigma*(1-(Se+Re)/Nemax)*Re - d*((C/(C+Cr50))*Re);
    dSm = lambda*(1-(Sm+Rm)/Nmmax)*Sm -eta*((L/(L+Ls50))*Sm)-phi*Sm;
    dRm = (1-p)*lambda*(1-(Sm+Rm)/Nmmax)*Rm - eta*((L/(L+Lr50))*Rm)-phi*Rm;
    dY=c(dS, dP, dC, dL, dSe,dRe,dSm,dRm);
    list(dY, signal= dose)
  }
  )
}

parms=c(k=0.44, beta= 0.0016, alpha= 0.16, vartheta= 0.53, sigma=0.28, Nemax=484577,
        Cs50= 0.25, Cr50= 0.471, d= 1.01, gamma=0.16 , lambda= 0.25, Nmmax=512861 ,
        Ls50=0.016 , Lr50= 0.125,delta= 0.24, eta=0.40, c=0.05, p=0.05, phi=0.14)
xstart= c(S=0,P=0,C=0,L=0,Se=400000,Re=1000,Sm=40000,Rm=1000)

Single_high1 = as.data.frame ( lsoda(xstart, times, optimalmodel, parms))




combined_data <- rbind(Single_high1, Multiple_low, Single_high)
combined_data$dataset <- rep(c("7.5mg/kg_single", "5mg/kg_three", "12.5mg/kg_single"), each = nrow(Single_high1))

# Create the line plot
plot1<-ggplot(combined_data, aes(x = time, y = Se, color = dataset)) +
  geom_line() + theme_minimal()+  
  labs(title = ("(a)"),
       x = "Time(hr)",
       y = expression("Susceptible"~italic("E. coli"))) +
  scale_color_manual(values = c("red", "blue", "green")) +
  scale_y_continuous(labels = scales::comma_format(scale = 1))+
  theme (panel.grid=element_blank(),panel.border = element_rect(colour = "gray", fill = NA, size = 1), plot.title=element_text(face="bold", size=12, vjust=2, family="Arial"))


# Create the line plot
plot2<-ggplot(combined_data, aes(x = time, y = Re, color = dataset)) +
  geom_line() + theme_minimal()+  
  labs(title = ("(b)"), 
       x = "Time(hr)",
       y = expression("Resistant"~italic("E. coli"))) +
  scale_color_manual(values = c("red", "blue", "green")) +
  scale_y_continuous(labels = scales::comma_format(scale = 1))+
  theme (panel.grid=element_blank(),panel.border = element_rect(colour = "gray", fill = NA, size = 1), plot.title=element_text(face="bold", size=12, vjust=2, family="Arial"))


combined_plot1 <- plot1 + plot2 +   # we used library patchwork to combine plot into one plot
  plot_layout(nrow = 2, heights = c(1, 1))  # Adjust layout if needed

# Print the combined plot
print(combined_plot1)

#Saving plot with specific width, height and high resolution
ggsave("E.coli.png", plot = combined_plot1, dpi = 300, width = 5, height = 4)

# Create the line plot
plot3<-ggplot(combined_data, aes(x = time, y = Sm, color = dataset)) +
  geom_line() + theme_minimal() + # to make shadow background then just remove the theme minimal part
  labs(title = ("(a)"),
       x = "Time(hr)",
       y = expression("Susceptible"~italic("P. multocida"))) +
  scale_color_manual(values = c("red", "blue", "green"))+
  scale_y_continuous(labels = scales::comma_format(scale = 1))+
  theme (panel.grid=element_blank(),panel.border = element_rect(colour = "gray", fill = NA, size = 1), plot.title=element_text(face="bold", size=12, vjust=2, family="Arial"))


# Create the line plot
plot4<-ggplot(combined_data, aes(x = time, y = Rm, color = dataset)) +
  geom_line() + theme_minimal() + # to make shadow background then just remove the theme minimal part
  labs(title = ("(b)"),
       x = "Time(hr)",
       y = expression("Resistant"~italic("P. multocida"))) +
  scale_color_manual(values = c("red", "blue", "green")) +
  scale_y_continuous(labels = scales::comma_format(scale = 1))+
  theme (panel.grid=element_blank(),panel.border = element_rect(colour = "gray", fill = NA, size = 1), plot.title=element_text(face="bold", size=12, vjust=2, family="Arial"))


combined_plot2 <- plot3 + plot4 +   # we used library patchwork to combine plot into one plot
  plot_layout(nrow = 2, heights = c(1, 1))  # Adjust layout if needed

# Print the combined plot
print(combined_plot2)

#Saving plot with specific width, height and high resolution
ggsave("P.multocida.png", plot = combined_plot2, dpi = 300, width = 5, height = 4)



# Model with proposed scenarios
# 6.25mg/kg- twice

times <- seq(0, 150, by = 1)
signa <- data.frame(times = times, import = rep(0, length(times)))

signa$import[signa$times==1]<-1250
signa$import[signa$times==25]<-1250
input <- approxfun(signa, rule = 2)
input(c(0,1,1.95833333,2,2.04166667,3,24,125))

optimalmodel = function(t,x,parms){
  with(as.list(c(parms,x)), {
    
    dose= input(t) 
    
    dS = -k*S+dose;
    dP = beta*k*S - (alpha+gamma)*P;
    dC = alpha*P - vartheta*C;
    dL = gamma*P - delta*L;
    # Check conditions and set Sm and Rm to 0 if less than 1
    if (Sm < 1) Sm = 0
    if (Rm < 1) Rm = 0
    dSe = sigma*(1-(Se+Re)/Nemax)*Se -d*((C/(C+Cs50))*Se);
    dRe = (1-c)*sigma*(1-(Se+Re)/Nemax)*Re - d*((C/(C+Cr50))*Re);
    dSm = lambda*(1-(Sm+Rm)/Nmmax)*Sm -eta*((L/(L+Ls50))*Sm)-phi*Sm;
    dRm = (1-p)*lambda*(1-(Sm+Rm)/Nmmax)*Rm - eta*((L/(L+Lr50))*Rm)-phi*Rm;
    dY=c(dS, dP, dC, dL, dSe,dRe,dSm,dRm);
    list(dY, signal= dose)
  }
  )
}

parms=c(k=0.44, beta= 0.0016, alpha= 0.16, vartheta= 0.53, sigma=0.28, Nemax=484577,
        Cs50= 0.25, Cr50= 0.471, d= 1.01, gamma=0.16 , lambda= 0.25, Nmmax=512861 ,
        Ls50=0.016 , Lr50= 0.125,delta= 0.24, eta=0.40, c=0.05, p=0.05, phi=0.14)
xstart= c(S=0,P=0,C=0,L=0,Se=100000,Re=1000,Sm=100000,Rm=1000)

Propose1 = as.data.frame ( lsoda(xstart, times, optimalmodel, parms))



# 7.5 mg/kg-twice

times <- seq(0, 150, by = 1)
signa <- data.frame(times = times, import = rep(0, length(times)))

signa$import[signa$times==1]<-1500
signa$import[signa$times==25]<-1500
input <- approxfun(signa, rule = 2)
input(c(0,1,1.95833333,2,2.04166667,3,24,125))

optimalmodel = function(t,x,parms){
  with(as.list(c(parms,x)), {
    
    dose= input(t) 
    
    dS = -k*S+dose;
    dP = beta*k*S - (alpha+gamma)*P;
    dC = alpha*P - vartheta*C;
    dL = gamma*P - delta*L;
    # Check conditions and set Sm and Rm to 0 if less than 1
    if (Sm < 1) Sm = 0
    if (Rm < 1) Rm = 0
    dSe = sigma*(1-(Se+Re)/Nemax)*Se -d*((C/(C+Cs50))*Se);
    dRe = (1-c)*sigma*(1-(Se+Re)/Nemax)*Re - d*((C/(C+Cr50))*Re);
    dSm = lambda*(1-(Sm+Rm)/Nmmax)*Sm -eta*((L/(L+Ls50))*Sm)-phi*Sm;
    dRm = (1-p)*lambda*(1-(Sm+Rm)/Nmmax)*Rm - eta*((L/(L+Lr50))*Rm)-phi*Rm;
    dY=c(dS, dP, dC, dL, dSe,dRe,dSm,dRm);
    list(dY, signal= dose)
  }
  )
}

parms=c(k=0.44, beta= 0.0016, alpha= 0.16, vartheta= 0.53, sigma=0.28, Nemax=484577,
        Cs50= 0.25, Cr50= 0.471, d= 1.01, gamma=0.16 , lambda= 0.25, Nmmax=512861 ,
        Ls50=0.016 , Lr50= 0.125,delta= 0.24, eta=0.40, c=0.05, p=0.05, phi=0.14)
xstart= c(S=0,P=0,C=0,L=0,Se=100000,Re=1000,Sm=100000,Rm=1000)

Propose2 = as.data.frame ( lsoda(xstart, times, optimalmodel, parms))
#out =  ode(y = xstart, times = times, func = optimalmodel, parms)

# 3.75 mg/kg-twice

times <- seq(0, 150, by = 1)
signa <- data.frame(times = times, import = rep(0, length(times)))

signa$import[signa$times==1]<-750
signa$import[signa$times==25]<-750
input <- approxfun(signa, rule = 2)
input(c(0,1,1.95833333,2,2.04166667,3,24,125))

optimalmodel = function(t,x,parms){
  with(as.list(c(parms,x)), {
    
    dose= input(t) 
    
    dS = -k*S+dose;
    dP = beta*k*S - (alpha+gamma)*P;
    dC = alpha*P - vartheta*C;
    dL = gamma*P - delta*L;
    # Check conditions and set Sm and Rm to 0 if less than 1
    if (Sm < 1) Sm = 0
    if (Rm < 1) Rm = 0
    dSe = sigma*(1-(Se+Re)/Nemax)*Se -d*((C/(C+Cs50))*Se);
    dRe = (1-c)*sigma*(1-(Se+Re)/Nemax)*Re - d*((C/(C+Cr50))*Re);
    dSm = lambda*(1-(Sm+Rm)/Nmmax)*Sm -eta*((L/(L+Ls50))*Sm)-phi*Sm;
    dRm = (1-p)*lambda*(1-(Sm+Rm)/Nmmax)*Rm - eta*((L/(L+Lr50))*Rm)-phi*Rm;
    dY=c(dS, dP, dC, dL, dSe,dRe,dSm,dRm);
    list(dY, signal= dose)
  }
  )
}

parms=c(k=0.44, beta= 0.0016, alpha= 0.16, vartheta= 0.53, sigma=0.28, Nemax=484577,
        Cs50= 0.25, Cr50= 0.471, d= 1.01, gamma=0.16 , lambda= 0.25, Nmmax=512861 ,
        Ls50=0.016 , Lr50= 0.125,delta= 0.24, eta=0.40, c=0.05, p=0.05, phi=0.14)
xstart= c(S=0,P=0,C=0,L=0,Se=100000,Re=1000,Sm=100000,Rm=1000)

Propose3 = as.data.frame ( lsoda(xstart, times, optimalmodel, parms))
#out =  ode(y = xstart, times = times, func = optimalmodel, parms)


# 12.5 mg/kg-twice

times <- seq(0, 150, by = 1)
signa <- data.frame(times = times, import = rep(0, length(times)))

signa$import[signa$times==1]<-2500
signa$import[signa$times==25]<-2500
input <- approxfun(signa, rule = 2)
input(c(0,1,1.95833333,2,2.04166667,3,24,125))

optimalmodel = function(t,x,parms){
  with(as.list(c(parms,x)), {
    
    dose= input(t) 
    
    dS = -k*S+dose;
    dP = beta*k*S - (alpha+gamma)*P;
    dC = alpha*P - vartheta*C;
    dL = gamma*P - delta*L;
    # Check conditions and set Sm and Rm to 0 if less than 1
    if (Sm < 1) Sm = 0
    if (Rm < 1) Rm = 0
    dSe = sigma*(1-(Se+Re)/Nemax)*Se -d*((C/(C+Cs50))*Se);
    dRe = (1-c)*sigma*(1-(Se+Re)/Nemax)*Re - d*((C/(C+Cr50))*Re);
    dSm = lambda*(1-(Sm+Rm)/Nmmax)*Sm -eta*((L/(L+Ls50))*Sm)-phi*Sm;
    dRm = (1-p)*lambda*(1-(Sm+Rm)/Nmmax)*Rm - eta*((L/(L+Lr50))*Rm)-phi*Rm;
    dY=c(dS, dP, dC, dL, dSe,dRe,dSm,dRm);
    list(dY, signal= dose)
  }
  )
}

parms=c(k=0.44, beta= 0.0016, alpha= 0.16, vartheta= 0.53, sigma=0.28, Nemax=484577,
        Cs50= 0.25, Cr50= 0.471, d= 1.01, gamma=0.16 , lambda= 0.25, Nmmax=512861 ,
        Ls50=0.016 , Lr50= 0.125,delta= 0.24, eta=0.40, c=0.05, p=0.05, phi=0.14)
xstart= c(S=0,P=0,C=0,L=0,Se=100000,Re=1000,Sm=100000,Rm=1000)

Propose4 = as.data.frame ( lsoda(xstart, times, optimalmodel, parms))
#out =  ode(y = xstart, times = times, func = optimalmodel, parms)


# 4.15 mg/kg-three

times <- seq(0, 150, by = 1)
signa <- data.frame(times = times, import = rep(0, length(times)))

signa$import[signa$times==1]<-1666.66
signa$import[signa$times==25]<-1666.66
signa$import[signa$times==49]<-1666.66
input <- approxfun(signa, rule = 2)
input(c(0,1,1.95833333,2,2.04166667,3,24,125))

optimalmodel = function(t,x,parms){
  with(as.list(c(parms,x)), {
    
    dose= input(t) 
    
    dS = -k*S+dose;
    dP = beta*k*S - (alpha+gamma)*P;
    dC = alpha*P - vartheta*C;
    dL = gamma*P - delta*L;
    # Check conditions and set Sm and Rm to 0 if less than 1
    if (Sm < 1) Sm = 0
    if (Rm < 1) Rm = 0
    dSe = sigma*(1-(Se+Re)/Nemax)*Se -d*((C/(C+Cs50))*Se);
    dRe = (1-c)*sigma*(1-(Se+Re)/Nemax)*Re - d*((C/(C+Cr50))*Re);
    dSm = lambda*(1-(Sm+Rm)/Nmmax)*Sm -eta*((L/(L+Ls50))*Sm)-phi*Sm;
    dRm = (1-p)*lambda*(1-(Sm+Rm)/Nmmax)*Rm - eta*((L/(L+Lr50))*Rm)-phi*Rm;
    dY=c(dS, dP, dC, dL, dSe,dRe,dSm,dRm);
    list(dY, signal= dose)
  }
  )
}

parms=c(k=0.44, beta= 0.0016, alpha= 0.16, vartheta= 0.53, sigma=0.28, Nemax=484577,
        Cs50= 0.25, Cr50= 0.471, d= 1.01, gamma=0.16 , lambda= 0.25, Nmmax=512861 ,
        Ls50=0.016 , Lr50= 0.125,delta= 0.24, eta=0.40, c=0.05, p=0.05, phi=0.14)
xstart= c(S=0,P=0,C=0,L=0,Se=100000,Re=1000,Sm=100000,Rm=1000)

Propose5 = as.data.frame ( lsoda(xstart, times, optimalmodel, parms))
#out =  ode(y = xstart, times = times, func = optimalmodel, parms)



combined_data1 <- rbind(Propose1, Propose2, Propose3, Propose4, Propose5)
combined_data1$Scenario <- rep(c("6.25mg/kg_twice", "7.5mg/kg_twice", "3.75mg/kg_twice","12.5mg/kg_twice","4.15mg/kg_three"), each = nrow(Propose1))


# Create the line plot
plot5<-ggplot(combined_data1, aes(x = time, y = Se, color = Scenario)) +
  geom_line() + theme_minimal() + # to make shadow background then just remove the theme minimal part
  labs(title = ("(a)"),
       x = "Time(hr)",
       y = expression(italic("Se"))) +  # Used expression function here to make one work italic in the y axis legend  
  scale_color_manual(values = c("red", "blue", "green","orange","Brown")) +
  scale_y_continuous(labels = scales::comma_format(scale = 1))+
  theme (panel.grid=element_blank(),panel.border = element_rect(colour = "gray", fill = NA, size = 1), 
         plot.title=element_text(face="bold", size=12, vjust=2, family="Arial"),
         legend.position = "none")


# Create the line plot
plot6<-ggplot(combined_data1, aes(x = time, y = Re, color = Scenario)) +
  geom_line() + theme_minimal() +
  labs(title = ("(b)"),
       x = "Time(hr)",
       y = expression(italic("Re"))) +
  scale_color_manual(values = c("red", "blue", "green","orange","Brown")) +
  scale_y_continuous(labels = scales::comma_format(scale = 1))+
  theme (panel.grid=element_blank(),panel.border = element_rect(colour = "gray", fill = NA, size = 1), 
         plot.title=element_text(face="bold", size=12, vjust=2, family="Arial"),
         legend.position = "none")

combined_plot3 <- plot5 + plot6 +   # we used library patchwork to combine plot into one plot
  plot_layout(nrow = 1, heights = c(1, 1))  # Adjust layout if needed

# Print the combined plot
print(combined_plot3)

#Saving plot with specific width, height and high resolution
ggsave("final.E.coli.proposed.png", plot = combined_plot3, dpi = 300, width = 8, height = 4)



# Create the line plot
plot7<-ggplot(combined_data1, aes(x = time, y = Sm, color = Scenario)) +
  geom_line() + theme_minimal() +
  labs(title = ("(c)"),     
       x = "Time(hr)",
       y =  expression(italic("Sm"))) +
  scale_color_manual(values = c("red", "blue", "green","orange","Brown")) +
  scale_y_continuous(labels = scales::comma_format(scale = 1))+
  theme (panel.grid=element_blank(),panel.border = element_rect(colour = "gray", fill = NA, size = 1), 
         plot.title=element_text(face="bold", size=12, vjust=2, family="Arial"),
         legend.position = "none")



# Create the line plot
plot8<-ggplot(combined_data1, aes(x = time, y = Rm, color = Scenario)) +
  geom_line() + theme_minimal() + # to make shadow background then just remove the theme minimal part
  labs(title = ("(d)"),    
       x = "Time(hr)",
       y =  expression(italic("Rm"))) +
  scale_color_manual(values = c("red", "blue", "green","orange","Brown")) +
  scale_y_continuous(labels = scales::comma_format(scale = 1))+
  theme (panel.grid=element_blank(),panel.border = element_rect(colour = "gray", fill = NA, size = 1), 
         plot.title=element_text(face="bold", size=12, vjust=2, family="Arial"),
         legend.position = "bottom", 
         legend.direction = "horizontal",
         legend.justification = "right")

combined_plot4 <- (plot5+plot7)/(plot6+ plot8) +   # we used library patchwork to combine plot into one plot
  plot_layout(nrow = 2, heights = c(1, 1))  # Adjust layout if needed

# Print the combined plot
print(combined_plot4)

#Saving plot with specific width, height and high resolution
ggsave("P.multocida.proposed.png", plot = combined_plot4, dpi = 300, width = 6, height = 6)
