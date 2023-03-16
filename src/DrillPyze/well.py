import well_profile as wp
import numpy as np

class wellbore:
    """
    A class for a wellbore, which includes a well path, a casing
    depth and a dimension.  This class uses the well_profile library to 
    interpolate along a wellbore survey.

    """


    def __init__(self, well_survey, casing_depth, hole_size):
        """
        Initialize a new Wellbore object.
        
        Parameters
        ----------
        survey : dict
            A dictionary containing one of the following keys:value pairs:
            ('xlsx' or 'csv':'path_to_survey') or ('dict':{'MD':numpy array of MD, 'INC':numpy array of INC, 'AZ':numpy array of Az})
        casing_depth : float
            Casing shoe depth
        hole_size : float
            hole size
            
        """

        if type(well_survey) is dict:
            if 'xlsx' in well_survey:
                self.well_trajectory = wp.load(well_survey['xlsx'], set_info={'units':'imperial'}, set_start={'north':0, 'east':0, 'depth':0})
            elif 'csv' in well_survey:
                self.well_trajectory = wp.load(well_survey['csv'], set_info={'units':'imperial'}, set_start={'north':0, 'east':0, 'depth':0})
            else:
                print('Not Implemented')
                self.well_trajectory = None

        self.casing_depth = casing_depth
        self.hole_size = hole_size

    def get_trajectory(self):
        """
        Return the well trajectory object
        """
        return self.well_trajectory

    def get_interpolated(self, md):
        """
        Interpolated wellpath for the depths defined.

        Parameters
        ==========
        md : numpy array
            Array of interpolation points
        """
        if self.well_trajectory is None:
            return -1
        else:
            tvd = np.zeros(len(md))
            inc = np.zeros(len(md))
            azi = np.zeros(len(md))

            # Iterate through the MDs to get the interpolated well trajectory at each point
            for i in range(0,len(md)):
                pt = self.well_trajectory.get_point(md[i])
                tvd[i] = pt['tvd']

                # well_profile does not define 'inc' and 'azi' consistently
                try:
                    inc[i] = pt['inc']
                except:
                    if i > 0:
                        inc[i] = inc[i-1]
                    else:
                        inc[i] = 0
                try:
                    
                    azi[i] = pt['azi']
                except:
                    if i > 0:
                        azi[i] = azi[i-1]
                    else:
                        azi[i] = 0  

            self.md = md
            self.tvd = tvd
            self.inc = inc
            self.azi = azi