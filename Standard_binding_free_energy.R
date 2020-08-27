library(zoo)

########################################################################################################################
########################################################################################################################
# caculate the binding free energy
#
# All your PMFs should be specified in the PMFs section below
# The results will be written in the file: binding_free_energy.nrj
#Author: Emmanuel Edouard Moutoussamy
#
#email: emmanuel.moutoussamy@uib.no
#
# 2nd of december 2017
#
# If you want mor information about the method, please see: https://pubs.acs.org/doi/10.1021/ct3008099
#
#
########################################################################################################################
########################################################################################################################

outfile = "binding_free_energy.nrj"
write("#Contribution numerator(N) denominator(D) Ratio(N/D) Energy_(kcal/mol)",file = outfile) #write header

########################################################################################################################
# Parameters: modified it according to your simulations parameters
########################################################################################################################

kb = 0.0019872041 # in kcal/K/mol
beta = 1.677399
T = 310 # Kelvin
k_conf = 15 #force constant used for RMSD
k_ori = 0.1 #force constant used for restraining the orientation
k_pos = 0.1 #force constant used for restraining the position
r0_conf_bound = 0 #equilibrium value for the prot/lig/ect conf
r0_conf_unbound = 0 #equilibrium value for the prot/lig/ect conf
r0_ori_theta = 10 #Equilibrium value for the Theta angle
r0_ori_Phi = 16 #Equilibrium value for the Phi angle
r0_ori_Psi = 15 #Equilibrium value for the Psi angle
r0_pol_theta = 29 #Equilibrium value for the theta polar angle
r0_pol_phi = -39 #Equilibrium value for the phi polar angle
r_star = 23 #defined r* (distance where to two partnes are not interacting anymore)

########################################################################################################################
# PMFs: File the PMFs for each contributions
########################################################################################################################


bound = 'Rmsdbound0.czar.pmf'
Theta_euler = 'merge_Theta.grad.pmf'
Phi_euler = 'merge_Phi.grad.pmf'
Psi_euler = 'merge_psi.grad.pmf'
theta_polar = 'merge_theta2.grad.pmf'
phi_polar = 'PolarPhi_hs_use.pmf'
sep = 'merge_05_03_20.pmf'
unbound = 'RmsdUnbound0.czar.pmf'

########################################################################################################################
# LETS DO THE CALCULATION NOW :)
########################################################################################################################

########################################################################################################################
# Function for contribution at site
########################################################################################################################

contrib <- function(k,r0,pmf){
        #Function to calculate the contribution of restraining the mobile partner in the binding site
        # works for RMSD bound, Theta, Phi and Psi
        numerator = read.table(pmf) #read the data
        denominator = numerator #duplicate the data
        u = (k*(denominator[,1] - r0)^2)/2
        denominator[,2] = denominator[,2] + u
        numerator[,2] = exp(-beta*numerator[,2])
        denominator[,2] = exp(-beta*denominator[,2])
        x = seq(dim(denominator)[1])
        id <- order(x)
        AUC_denominator <- sum(diff(denominator[id,1])*rollmean(denominator[id,2],2)) #integration of the numerator
        AUC_numerator <- sum(diff(numerator[id,1])*rollmean(numerator[id,2],2))#integration of the denominator
        exp_beta = AUC_numerator/AUC_denominator
        G = log(exp_beta)/beta

        results = c(AUC_numerator,AUC_denominator,exp_beta,G)
        return(results) # return denominator,numerator, ratio and energy
}

########################################################################################################################
# Function for writting the results in outfile
########################################################################################################################


write_results <- function(contribution,numerator,denominator,ratio,nrj){
        # Write results in a csv format in the file "binding_free_energy.nrj"
        # Easy to parse with python: if "#" not line line , line = line.split()
        # Each line that not contain "#" has 6 elements (contrib, numerator, denominator, ratio and G)
        # The final binding free energy is indicated with the tag #GBIND
        # The Keq is indicated with the #KEQ tag
        output_line = paste(contribution,numerator,denominator,ratio,nrj,sep = " ")
        write(output_line,file=outfile,append = TRUE)
}


########################################################################################################################
# Contribution for restraining the conformation of the prot/lig/ect (mobile patner)
########################################################################################################################


ConfBound = contrib(k_conf,r0_conf_bound,bound)
G_RMSDb = ConfBound[4]
write_results("Rmsdb",ConfBound[1],ConfBound[2],ConfBound[3],ConfBound[4])

########################################################################################################################
# Contribution for restraining the orientation of the prot/lig/ect (mobile patner) with respect to Theta (Euler angle)
########################################################################################################################

Theta = contrib(k_ori,r0_ori_theta,Theta_euler)
G_Theta = Theta[4]
write_results("Theta",Theta[1],Theta[2],Theta[3],Theta[4])

########################################################################################################################
# Contribution for restraining the orientation of the prot/lig/ect (mobile patner) with respect to Phi (Euler angle)
########################################################################################################################

Phi = contrib(k_ori,r0_ori_Phi,Phi_euler)
G_phi = Phi[4]
write_results("Phi",Phi[1],Phi[2],Phi[3],Phi[4])

########################################################################################################################
# Contribution for restraining the orientation of the prot/lig/ect (mobile patner) with respect to Psi (Euler angle)
########################################################################################################################

Psi = contrib(k_ori,r0_ori_Psi,Psi_euler)
G_Psi = Psi[4]
write_results("Psi",Psi[1],Psi[2],Psi[3],Psi[4])

########################################################################################################################
# Contribution for restraining the position of the prot/lig/ect (mobile patner) with respect to theta (Polar angle)
########################################################################################################################

pol_theta = contrib(k_pos,r0_pol_theta,theta_polar)
G_theta = pol_theta[4]
write_results("theta",pol_theta[1],pol_theta[2],pol_theta[3],pol_theta[4])


########################################################################################################################
# Contribution for restraining the position of the prot/lig/ect (mobile patner) with respect to phi (Polar angle)
########################################################################################################################

pol_phi = contrib(k_pos,r0_pol_phi,phi_polar)
G_phi = pol_phi[4]
write_results("phi",pol_phi[1],pol_phi[2],pol_phi[3],pol_phi[4])


########################################################################################################################
# Contribution for releasing restraints on Theta, Phi and Psi (Euler angles) on bulk
########################################################################################################################


#Function regarding Theta
k_Theta = k_ori *(180/pi)^2
Theta0 = r0_Ori_theta * (pi/180)
Theta <- function(x){
        sin(x)*exp(-beta*(0.5*k_Theta)*((x-Theta0))^2)
}

#Function regarding Phi
k_Phi = k_ori *(180/pi)^2
Phi0 = r0_ori_Phi * (pi/180)
Phi <- function(x){
        exp(-beta*(0.5*k_Phi)*((x-Phi0))^2)
}


#Function regarding Psi
k_Psi = k_ori *(180/pi)^2
Psi0 = r0_ori_Psi * (pi/180)
Psi <- function(x){
        exp(-beta*(0.5*k_Psi)*((x-Psi0))^2)
}



# Integration of the three parts
Theta_part = integrate(Theta, lower= 0,upper =pi)$value
Phi_part = integrate(Phi, lower= 0,upper =2*pi)$value
Psi_part = integrate(Psi, lower= 0,upper =2*pi)$value


exp_beta_G_bulk_o = (1/(8*pi^2)) * Theta_part * Phi_part * Psi_part
G_bulk_o = log(exp_beta_G_bulk_o)/-beta


write_results("Theta/Phi/Psi_bulk","-","-","-",G_bulk_o)


########################################################################################################################
# I* term (Sep)
########################################################################################################################


sep = read.table(sep)

# /!\  remove Jacobian correction. comment the following line if you use 'HideJacobian on' during your calculation /!\
sep[,2] = sep[,2]  + (2*kb*T*log10(sep[,1]))



r_star_nearest = which.min(abs(sep[,1]-r_star))

sep[,2] = sep[,2] - sep[r_star_nearest,2]
sep[,2] = exp(-beta*sep[,2])

x = seq(dim(sep)[1])
id <- order(x)

AUC_I_star <- sum(diff(sep[id,1])*rollmean(sep[id,2],2))


write_results("I*","-","-","-",AUC_I_star)


########################################################################################################################
# S*  term
########################################################################################################################

# effectively describes the small fraction of surface area on the sphere of radius r* centered on the binding site that
# is accessible to the restrained protein/lig/ect (mobile partner)


# Polar theta part
k_theta = k_pos *(180/pi)^2
theta0 = r0_pol_theta * (pi/180)
theta <- function(x){
  sin(x)*exp(-beta*(0.5*k_theta)*((x-theta0))^2)
}

#Polar phi part
k_phi = k_pos *(180/pi)^2
phi0 = (r0_pol_phi+360) * (pi/180)
phi <- function(x){
  exp(-beta*(0.5*k_phi)*((x-phi0))^2)
}


theta_part = integrate(theta, lower= 0,upper =pi/2)$value
phi_part = integrate(phi, lower= 0,upper =2*pi)$value

S_star = ((r_star)^2) * theta_part * phi_part

write_results("S*","-","-","-",S_star)


########################################################################################################################
# Contribution for releasing restrains on the protein/lig/ect (mobile partner) conf. in bulk
########################################################################################################################



numerator = read.table(unbound)
denominator = numerator


u = (k_conf *(denominator[,1]- r0_conf_unbound)^2)/2


denominator[,2] = numerator[,2] + u
numerator[,2] = exp(-beta*numerator[,2])
denominator[,2] = exp(-beta*denominator[,2])


x = seq(dim(denominator)[1])
id <- order(x)

AUC_denominator <- sum(diff(denominator[id,1])*rollmean(denominator[id,2],2))
AUC_numerator <- sum(diff(numerator[id,1])*rollmean(numerator[id,2],2))

exp_beta_G_free = AUC_numerator/AUC_denominator
G_RMSDf = log(exp_beta_G_free)/beta


write_results("RMSDf",AUC_numerator,AUC_denominator,exp_beta_G_free,G_RMSDf)

########################################################################################################################
########################################################################################################################
# FINAL CALCULATION OF THE BINDING FREE ENERGY
########################################################################################################################
########################################################################################################################


SI_star = S_star*AUC_I_star #S*I* term
exp_contrib = exp(-beta*(G_RMSDf+G_bulk_o-(G_phi+G_theta)-(G_Theta+G_Phi+G_Psi)-G_RMSDb))

Keq = SI_star*exp_contrib #equilibium constant
output_line = paste("\n#KEQ",Keq)
write(output_line,file=outfile,append = TRUE)

G_bind = -kb*T*log(Keq*(1/1661)) #binding free energy calculation

output_line = paste("\n#GBIND ",G_bind)
write(output_line,file=outfile,append = TRUE)



################################# THAT IS ALL FOLKS ! ##################################################################
########################################################################################################################