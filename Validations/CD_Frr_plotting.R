# Plot CD measurements for Frr
# Author: Monika Pepelnjak

# Goal of the script:
# Plot the scaled and fitted values for CD measurements

# Fitting and scaling is performed by python script

### CD measurements ### 
# Load in the raw values, the fitted data and the melting temperatures
CD_Frr <- read.csv("/Users/moni/Documents/Phd/Collaborations/C6/CD/Frr_scaled_y_values.csv")
CD_Fit <- read.csv("/Users/moni/Documents/Phd/Collaborations/C6/CD/Frr_Fitted_values.csv") %>%
  `colnames<-`(c("Temperature", "Control", "TMAO"))
CD_Melt <- read.csv("/Users/moni/Documents/Phd/Collaborations/C6/CD/Frr_output_parameters.csv") %>%
  mutate(Tm = Tm-273.15) 

# Reshape in the right format for ggplot
CD_Frr %<>%
  reshape2::melt(., id.var="Temperature") %>%
  `colnames<-`(c("Temperature","Condition","Intensity" ))
CD_Fit %<>%
  reshape2::melt(., id.var="Temperature") %>%
  `colnames<-`(c("Temperature", "Condition", "Intensity"))

# Plot the data
ggplot() +
  geom_point(data = CD_Frr, aes(x=Temperature-273.15, y=Intensity, col=Condition)) + 
  geom_line(data = CD_Fit, aes(x=Temperature-273.15, y=Intensity, col=Condition)) + 
  scale_color_manual(values=c("Control" = "azure3", my_palette["TMAO"])) +
  theme_classic() +
  ylab("Fraction unfolded") +
  xlab("Temperature")
  

