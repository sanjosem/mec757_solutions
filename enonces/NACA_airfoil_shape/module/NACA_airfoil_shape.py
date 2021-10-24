import numpy as np

class NACA_Airfoil:
    """
    Contain NACA airfoil parameters and associated functions to compute the form
    """
    def __init__(self, number, chord = 1.0,sharp=True):
        """
        Sets the airfoil parameter.
        
        Parameters
        ----------
        number: 4 digit string 
            NACA number
        chord : float
            chord of the airfoil
        sharp: boolean
            sharp or blunt TE (default sharp=True)

        """
        self.number = number
        self.camber_max = int(number[0])/100.
        self.camber_pos = int(number[1])/10.
        self.thickness = int(number[2:4])/100.
        self.chord = chord
        self.coeff = np.array((0.2966,-0.1260,-0.3516,0.2843,-0.1015))
        self.expo = np.array((0.5,1,2,3,4))
        if sharp:
            self.coeff[-1]=-0.1033

    def half_thickness(self,x):
        """
        Compute the half thickness of symmetric NACA 
        
        Parameters
        ----------
        x: coordinates along the chord
        """        
        xc = x / self.chord
        def yt(xc):
            y_t = np.zeros_like(xc)
            for ic in range(5):
                y_t +=  self.coeff[ic] * xc**(self.expo[ic])
            y_t *= self.thickness / 0.2
            return y_t

        self.shape = yt(xc) * self.chord


    def camber(self,x):
        """
        Compute the mean camber line of NACA 4 digits
        
        Parameters
        ----------
        x: coordinates along the chord
        """        
        xc = x / self.chord

        def yc(xc):
            y_c = np.zeros_like(xc)
            xc_up = xc < self.camber_pos
            xc_dn = xc >= self.camber_pos
            if np.any(xc_up):
                y_c[xc_up] =  self.camber_max / self.camber_pos**2 * (2.*self.camber_pos*xc[xc_up]-xc[xc_up]**2)
            y_c[xc_dn] =  self.camber_max / (1.-self.camber_pos)**2 * (1.-2.*self.camber_pos + 2.*self.camber_pos*xc[xc_dn]-xc[xc_dn]**2)
            
            return y_c

        self.mean_camber_line = yc(xc) * self.chord

    def deriv_camber(self,x):
        """
        Compute the mean camber line derivative of NACA 4 digits
        
        Parameters
        ----------
        x: coordinates along the chord
        """        
        xc = x / self.chord

        def dycdx(xc):
            dy_c = np.zeros_like(xc)
            xc_up = xc < self.camber_pos
            xc_dn = xc >= self.camber_pos
            if np.any(xc_up):
                dy_c[xc_up] =  2. * self.camber_max / self.camber_pos**2 * (self.camber_pos-xc[xc_up])
            dy_c[xc_dn] =  2.*self.camber_max / (1.-self.camber_pos)**2 * (self.camber_pos -xc[xc_dn])
            
            return dy_c

        self.deriv_camber_fct = dycdx(xc)

    def profile_side_shape(self,x):
        """
        Compute the suction side of NACA 4 digits
        
        Parameters
        ----------
        x: coordinates along the chord

        Returns
        -------
        xU: coordinates along the x-axis of the upper edge of the profile
        yU: coordinates along the y-axis of the upper edge of the profile
        xL: coordinates along the x-axis of the lower edge of the profile
        yL: coordinates along the y-axis of the lower edge of the profile
        """
        self.half_thickness(x)
        if self.camber_max>0:
            self.camber(x)
            self.deriv_camber(x)
            theta = np.arctan(self.deriv_camber_fct)
            xU = x - self.shape * np.sin(theta)
            yU = self.mean_camber_line + self.shape * np.cos(theta)
            xL = x + self.shape * np.sin(theta)
            yL = self.mean_camber_line - self.shape * np.cos(theta)
        else:
            xU = x
            yU = self.shape
            xL = x
            yL = - self.shape
        
        return xU,yU,xL,yL

