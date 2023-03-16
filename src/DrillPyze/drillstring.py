import pandas as pd
import numpy as np


class BHA:
    """
    A class describing a bottom hole assembly with the ability
    to discretize the string
    """

    def __init__(self, name, filename):
        """
        Create a new BHA based on the input excel file.  Note, the first component
        should be a single joint of drillpipe that will go to surface.

        Parameters
        ==========
        name : str
            Name of the BHA
        filename : str
            filename of the excel file containing the BHA.  Required columns are 'Component', 'Length', 'OD', 'ID', and 'Weight'.
            'Component' must be one of ('DP', 'HWDP', 'DC', 'Stab', 'Tool', 'Bit')
            'Weight' should be in lb/ft and is optional
        """
        self.name = name
        try:
            self.data = pd.read_excel(filename)

        except:
            print('Error loading BHA')


    def discretize(self, dx, bit_depth):
        """
        Discretize the drillstring into compoments of maximum length dx.  Not all components will be the same length
        
        Parameters
        ==========
        dx : float
            Maximal spacial discretization length
        bit_depth : float
            Desired length of the drillstring
        
        """
        self.bit_depth = bit_depth
        self.dx = dx

        l_dp = bit_depth - self.data.iloc[1:,1].sum()
        self.l_dp = l_dp
        n_elem = np.zeros(len(self.data), dtype=int)

        # Determine the number of discretization elements for each BHA component
        for i in range(0,len(self.data)):
            if i == 0:
                n_elem[i] = int(np.ceil(l_dp/dx))
            else:
                n_elem[i] = int(np.ceil(self.data.iloc[i,1]/dx))
        
        # Allocate numpy arrays
        elem = n_elem.sum()
        OD = np.zeros(elem)
        ID = np.zeros(elem)
        l = np.zeros(elem)
        x = np.zeros(elem)


        # Populate numpy arrays with OD, ID, J, dx and w
        offset = 0
        curr_length = 0
        for i in range(0,len(n_elem)):
            for j in range(0,n_elem[i]):
                OD[j+offset] = self.data.iloc[i,2]
                ID[j+offset] = self.data.iloc[i,3]
                if i==0:
                    l[j+offset] = l_dp/n_elem[i]
                else:
                    l[j+offset] = self.data.iloc[i,1]/n_elem[i]
                curr_length += l[j+offset]
                x[j+offset] = curr_length
            offset += n_elem[i]

        self.OD = OD
        self.ID = ID
        self.J = np.pi / 32 * (np.power(self.OD/12, 4) - np.power(self.ID/12, 4))
        self.l = l
        self.w = l * 2.67 * (np.power(self.OD, 2) - np.power(self.ID, 2))
        self.x = x
