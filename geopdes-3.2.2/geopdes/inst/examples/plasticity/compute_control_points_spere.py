import numpy as np
from igakit.cad import circle
from igakit.cad import revolve
from igakit.cad import ruled
from igakit.plot import plt


# ------------------------------------------------------------------------
# CLASS
# ------------------------------------------------------------------------

# solve Poisson problem in 1D
# div( grad(v)) = f
# standard isogeometric collocation method

class sphere_eighth:

    def __init__(self, R1 =.5, R2 =1.):
    
        # DATA
        # internal radius
        self.R1 = R1
        # external radius
        self.R2 = R2
        
    
    # ---------------------------------------------------------------------
    
    def _create_surface(self, R):    
    
        c1 = circle(radius=R, angle=(0,np.pi/4.))
        c1.rotate(np.pi/2, axis=0)
        surf = revolve(c1, point=0, axis=2, angle=[0,np.pi/2.])
        
        # print(c1.points)
        # plt.plot(c1, color='b')
        # plt.plot(surf, color='g')
  
        
        return surf
    
    # ---------------------------------------------------------------------
    
    def _create_sphere(self):
    
        surf_int = self._create_surface(self.R1)
        surf_ext = self._create_surface(self.R2)
        
        volume = ruled(surf_int, surf_ext)
        
        return volume
        
    # ---------------------------------------------------------------------
    
    def _refinement(self, vol, dir):
    
        vol.elevate(dir,1)
        
        return vol
        
    def _write_file(self, CP):
        
        f = open("geo_eighth_sphere.txt", "w")
        f.write("# nurbs mesh v.2.1\n")
        f.write("#\n")
        f.write("# 9-May-2024\n")
        f.write("#\n")
        f.write(" 3 3 1 0 1\n")
        f.write("PATCH 1 \n")
        f.write("   2   2   2\n")
        f.write("   3   3   3\n")
        f.write("0.0000000   0.0000000   0.0000000   1.0000000   1.0000000   1.0000000 \n")
        f.write("0.0000000   0.0000000   0.0000000   1.0000000   1.0000000   1.0000000 \n")
        f.write("0.0000000   0.0000000   0.0000000   1.0000000   1.0000000   1.0000000 \n")
        
        for j in [0,1,2,3]:
            for i in range(0,len(CP[j,:])):            
                f.write(str(np.round(CP[j,i], decimals=7))+'  ')
            f.write('\n')

        f.write('SUBDOMAIN 1\n')
        f.write('1')
        f.close()
    
    
        return
        
    # ---------------------------------------------------------------------
        
    def output_data(self):
    
        volume = self._create_sphere()
        volume = self._refinement(volume, 2)
        
        # print('degree: ')
        # print(volume.degree)
        # print('knots: ')
        # print(volume.knots)        
        plt.plot(volume, color='r')
        plt.show()
    
        control_points= volume.control   
        control_points = control_points.transpose().reshape((4,27), order='F')
        
        self._write_file(control_points)
        
        return 
    
# ------------------------------------------------------------------------
# RUN
# ------------------------------------------------------------------------    
    
if __name__ == '__main__':
    geometry = sphere_eighth(R1=100, R2=200)
    geometry.output_data()

    
  
