#This file will contain the functions to describe magnetospheric boundaries and functions to describe the emissivity. 

#You wouldn't run this file directly. 

import numpy as np 

def shue_func(theta, phi, r0, ay, az):
        '''This is the 3D Shue model defined in Jorgensen et al. (2019)
        
        Parameters
        ----------
        theta (rad) and phi (rad)
        r0 - subsolar magnetopause distance
        ay - alpha y parameter
        az - alpha z parameter 

        Returns
        -------
        r - radial distance at the angles theta and phi 
        '''

        ry = r0*((2/(1+np.cos(theta)))**ay)
        rz = r0*((2/(1+np.cos(theta)))**az)

        r = (ry*rz)/(((rz*np.cos(phi))**2 + (ry*np.sin(phi))**2)**0.5)

        return r 
        
def lin_scaled_func(theta, phi, a, beta_c, c, dn, ds, theta_n, theta_s, r0_lin, p0=1, p1=1, p2=1, p3=1):
        '''This function will work out r using the lin model. 
        
        Parameters
        ----------
        theta (rad) - Shue coords.
        phi (rad) - Shue coords. 
        a, beta_c, c, dn, ds, theta_n, theta_s, r0_lin - Lin coefficients in model. 
        #dipole - dipole tilt angle (rad)
        #pd - dynamic pressure in nPa
        #pm - magnetic pressure in nPa 
        #bz - IMF bz component in nT 
        p - parameter scaling factors. 
            p0 scales r0
            p1 scales flaring parameter beta 
            p2 scales indentation parameter Q (cusp depth) 
            p3 scales d in indentation shape (cusp shape/width)
            '''

        # Get coefficients if for some reason, they have not already been calculated. 
        #if r0_lin is None: 
        #    get_lin_coeffs(dipole, pd, pm, bz)
        
        # Get phi-n and phi-s.
        phi_n = np.arccos((np.cos(theta)*np.cos(theta_n)) + (np.sin(theta)*np.sin(theta_n)*np.cos(phi-(np.pi/2.))))
        phi_s = np.arccos((np.cos(theta)*np.cos(theta_s)) + (np.sin(theta)*np.sin(theta_s)*np.cos(phi-(3*np.pi/2.))))

        # Get f. 
        f = (np.cos(theta/2) + a[5]*np.sin(2*theta)*(1-np.exp(-theta)))**(p1*(beta_c[0] + beta_c[1]*np.cos(phi) + beta_c[2]*np.sin(phi) + beta_c[3]*(np.sin(phi)**2)))

        # Get Q. 
        Q = p2*c*np.exp(p3*dn*(phi_n**a[21])) + p2*c*np.exp(p3*ds*(phi_s**a[21]))

        # Get r. 
        r = p0*r0_lin*f + Q

        return r 
  
def get_lin_coeffs(dipole, pd, pm, bz):
        '''This gets the value of r0 in the Lin et al. (2010) model, which is a constant value 
        that depends on solar wind parameters. All of these functions are independent of beta and gamma. 
        
        Parameters
        ----------
        dipole
        pd
        pm
        bz
        
        Returns
        -------
        All coefficients are attached to self. 
        '''

        # Get a coefficients first. 
        a = np.array([12.544, -0.194, 0.305, 0.0573, 2.178, 0.0571, -0.999, 16.473, 0.00152, 0.382, 0.0431, -0.00763, -0.210, 0.0405, -4.430, -0.636, -2.600, 0.832, -5.328, 1.103, -0.907, 1.450])
        #self.a = a

        # Get beta coefficients - renamed delta. 
        beta_c = np.array([a[6] + a[7]*((np.exp(a[8]*bz) - 1)/(np.exp(a[9]*bz) + 1)), a[10], a[11] + a[12]*dipole, a[13]])
         
        # Get cn and cs coefficients (equal). 
        c = a[14]*(pd+pm)**a[15]

        # Get d coefficients. 
        dn = (a[16] + a[17]*dipole + a[18]*dipole**2)
        ds = (a[16] - a[17]*dipole + a[18]*dipole**2)
        
        # Get theta-n and theta-s coefficients.
        theta_n = a[19] + a[20]*dipole
        theta_s = a[19] - a[20]*dipole

        # Get the unscaled subsolar magnetopause radius. 
        r0_lin = 12.544*((pd+pm)**-0.194)*(1 + 0.305*((np.exp(0.0573*bz) -1 )/(np.exp(2.178*bz) + 1)))
        
        return a, beta_c, c, dn, ds, theta_n, theta_s, r0_lin 
        
def get_model_func(current_model):
        '''This will select the correct function for the desired model. '''
        
        if current_model == "jorg":
            def jorg_func(r, theta, phi, mp, bs, A1, A2, B, alpha, beta, ay_mp, az_mp, ay_bs, az_bs):
               
                '''This is the model from the Jorgensen paper. 
        
                Parameters
                ----------
                r - 3D array of r values.
                theta - 3D array of theta values. 
                phi - 3D array of phi values. 
                mp - subsolar magnetopause distance parameter
                bs - subsolar bowshock distance parameter
                A1 - parameter
                A2 - parameter
                B - parameter
                alpha - parameter
                beta - parameter
                ay_mp - ay magnetopause flaring parameter
                az_mp - az magnetopause flaring parameter
                ay_bs - ay bowshock flaring parameter
                az_bs - az bowshock flaring parameter
                '''

                eta = np.zeros(r.shape)

                # Calculate the radii to the magnetopause and bowshock for all 
                # combinations of theta and phi. 
                rmp = shue_func(theta, phi, mp, ay_mp, az_mp)
                rbs = shue_func(theta, phi, bs, ay_bs, az_bs)

                # Get indices inside MP, between MP and BS, and outside BS. 
                r1 = np.where(r < rmp)
                r2 = np.where((r >= rmp) & (r < rbs))
                r3 = np.where(r >= rbs)

                # Now calculate eta in each region. 
                eta[r1] = 0.0
                eta[r2] = (A1 + B*((np.sin(theta[r2]))**8))*((r[r2]/10)**(-alpha-(beta*(np.sin(theta[r2]))**2)))
                eta[r3] = A2*((r[r3]/10)**(-3))
        
                return eta
            return jorg_func
             
        elif current_model == "cmem":
            def cmem_func(r, theta, phi, a, beta_c, c, dn, ds, theta_n, theta_s, r0_lin, p0, bs, A1, A2, B, alpha, beta, p1, p2, p3, ay_bs, az_bs):
                '''
                This is the CMEM model, which will use the lin model to work out 
                the magnetopause location instead of the shue model. 

                Parameters
                ----------
                r - 3D array of r values.
                theta - 3D array of theta values. 
                phi - 3D array of phi values. 
                a, beta_c, c, dn, ds, theta_n, theta_s, r0_lin - Lin coefficients in model. 
                p0 - scaling factor on the subsolar magnetopause parameter 
                bs - subsolar bowshock distance parameter
                A1 - parameter
                A2 - parameter
                B - parameter
                alpha - parameter
                beta - parameter
                p1 - scaling factor on magnetopause flaring parameter
                p2 - scaling parameter on magnetopause indentation parameter 
                ay_bs - ay bowshock flaring parameter
                az_bs - az bowshock flaring parameter
                '''
            
                eta = np.zeros(r.shape)

                # Calculate the radii to the magnetopause and bowshock for all 
                # combinations of theta and phi. 
                rmp = lin_scaled_func(theta, phi, a, beta_c, c, dn, ds, theta_n, theta_s, r0_lin, p0, p1, p2, p3)
                rbs = shue_func(theta, phi, bs, ay_bs, az_bs)

                # Get indices inside MP, between MP and BS, and outside BS. 
                r1 = np.where(r < rmp)
                r2 = np.where((r >= rmp) & (r < rbs))
                r3 = np.where(r >= rbs)

                # Now calculate eta in each region. 
                eta[r1] = 0.0
                eta[r2] = A1*(np.exp(-B*(theta[r2]/2.)**4))*((r[r2]/10)**(-alpha-(beta*(np.sin(theta[r2]))**2)))
                eta[r3] = A2*((r[r3]/10)**(-3))
                
                return eta
            return cmem_func
        
        else:
            raise ValueError("{} not a valid model. 'jorg' or 'cmem' only atm.".format(current_model))

