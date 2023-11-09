######################################################################
# Author : Lars Berger                                               #
# Date : 24.10.2020                                                  #
#                                                                    #
# A virtual detector to take the full 3D velocity space coverage of  #
# STEP into account                                                  #
######################################################################

import numpy as np
#from pylib.etCoord import rotate
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.cm as cm

plt.rcParams["xtick.direction"] = "in"
plt.rcParams["xtick.top"] = "True"
plt.rcParams["ytick.direction"] = "in"
plt.rcParams["ytick.right"] = "True"
plt.rcParams["ytick.major.size"] = 6
plt.rcParams["xtick.major.size"] = 6
plt.rcParams["ytick.minor.size"] = 4
plt.rcParams["xtick.minor.size"] = 4
plt.rcParams["xtick.major.width"] = 2.
plt.rcParams["ytick.major.width"] = 2.
plt.rcParams["xtick.minor.width"] = 1.5
plt.rcParams["ytick.minor.width"] = 1.5
plt.rcParams["axes.linewidth"] = 2.
plt.rcParams["legend.frameon"] = False


def create_FoV_frame():
    fig = plt.figure(figsize=(11,10))
    fig.subplots_adjust(wspace=0.,hspace=0.)
    ax = []
    ax.append(fig.add_subplot(1,2,1))
    for c in range(4,7):
        for r in range(5):
            k = c+r*6
            ax.append(fig.add_subplot(5,6,k,sharex=ax[0],sharey=ax[0]))
    ax[0].set_xticks(np.arange(20,50.1,10.))
    ax[0].set_xticks(np.arange(20,50.1,2.),minor=True)
    ax[0].set_yticks(np.arange(-30,30.1,10.))
    ax[0].set_yticks(np.arange(-30,30.1,2.),minor=True)
    ax[0].set_xticklabels(['-20','-30','-40','-50'])
    ax[0].set_ylim(-35,35)
    ax[0].set_xlim(15,55)
    ax[0].set_ylabel(r'$\theta$ [°]',fontsize = 20, labelpad = -1.)
    for i in range(1,16):
        ax[i].set_xticks(np.arange(20,50.1,10.))
        ax[i].set_xticks(np.arange(20,50.1,2.),minor=True)
        ax[i].set_yticks(np.arange(-30,30.1,10.))
        ax[i].set_yticks(np.arange(-30,30.1,2.),minor=True)
        ax[i].set_xticklabels(['-20','-30','-40','-50'])
        ax[i].text(0.1,0.1,str(i),transform=ax[i].transAxes,fontsize=15)
        if i in [1,6,11]:
            ax[i].xaxis.tick_top()
        if i in [2,3,4,7,8,9,12,13,14]:
            plt.setp(ax[i].get_xticklabels(), visible=False) 
        if i in [3,6,7,8,9,10]:
            plt.setp(ax[i].get_yticklabels(), visible=False)
        if i in [11,12,13,14,15]:
            ax[i].yaxis.tick_right()
        if i in [1,2,3,4,5]:
            plt.setp(ax[i].get_yticklabels(), visible=False)
    fig.text(0.5,0.05,r'$\phi$ [°]',fontsize = 20)
    return fig,ax



class VDSTEP_num(object):
    def __init__(self, npy = 51, angstep = 1., B = np.array([np.cos(4/18*np.pi),-np.sin(4/18.*np.pi),0.0]),mag = False):
        """
        Coordinatesystem :
        x -> Axis through center of STEP (i.e. Pixel 8) orientation from pinhole to detector array (0,0,0 = center of Pinhole)
        y -> from center to background pixel (i.e. parallel to T in RTN)
        Z -> completes the right handed triade (i.e. parallel to N in RTN)
        """
        self.phi_d = np.arange(-16,16.1,angstep) # angles in Instrument Coordinates
        self.theta_d = np.arange(-30.,30.1,angstep) # angles in Instrument Coordinates
                                               # theta is not in spherical coordinates but against the main plane, i.e. centre of bckgrd and pixel 3,8,13
        self.phi_r = self.phi_d/180.*np.pi
        self.theta_r = self.theta_d/180.*np.pi
        self.phi_offs_d = -35.
        self.theta_offs_d = 0.
        self.phi_offs_r = self.phi_offs_d/180.*np.pi
        self.theta_offs_r = self.theta_offs_d/180.*np.pi
        self.d = 15 # Distance between pinhole(front) and plane of the ssd-array in mm
        self.B = B
        self.npy = npy
        self.npz = 2*self.npy
        self.dim = self.npy*self.npz
        self.directions = np.zeros((self.dim,3,self.phi_r.shape[0],self.theta_r.shape[0]))
        self.vdir = np.zeros((3,self.phi_r.shape[0],self.theta_r.shape[0]))
        self.cosmu = np.zeros((self.phi_r.shape[0],self.theta_r.shape[0]))
        self.mu = np.zeros((self.phi_r.shape[0],self.theta_r.shape[0]))
        self.A_pinhole = 0.14 * 0.28  # cm^2
        self.pA_pinhole = np.zeros((self.phi_r.shape[0],self.theta_r.shape[0]))
        self.hitfrac = np.zeros((16,self.phi_r.shape[0],self.theta_r.shape[0]))
        self.geometry_factor = np.zeros((16,self.phi_r.shape[0],self.theta_r.shape[0]))
        self.scale = np.ones((self.phi_r.shape[0],self.theta_r.shape[0]))*np.cos(self.theta_r)[np.newaxis,:]
        self.pixel_pos = np.zeros((16,3)) 
        self.mag = mag
        self._set_pixel_pos()
        self._set_pinhole()
        self._calc_pA()
        self._calc_directions()
        self._calc_geometry_factor()
        self._calc_vdir()
        self._calc_pitch_ang()
        #self._calc_pitch_ang_cov()
        
    def _calc_vdir(self):
        # calc unity velocity vectors in SUN-SC Coordinates
        for i,p in enumerate(self.phi_r):
            for j,t in enumerate(self.theta_r):
                th = -t + np.pi/2. + self.theta_offs_r
                ph = p + self.phi_offs_r
                self.vdir[0,i,j] = np.sin(th)*np.cos(ph)
                self.vdir[1,i,j] = np.sin(th)*np.sin(ph)
                self.vdir[2,i,j] = np.cos(th)

    def _calc_pitch_ang(self):
        Bx = self.B[0]/np.sqrt(self.B[0]**2+self.B[1]**2+self.B[2]**2)
        By = self.B[1]/np.sqrt(self.B[0]**2+self.B[1]**2+self.B[2]**2)
        Bz = self.B[2]/np.sqrt(self.B[0]**2+self.B[1]**2+self.B[2]**2)
        vx = self.vdir[0]
        vy = self.vdir[1]
        vz = self.vdir[2]
        vdim = vx.shape
        Bx = np.array([Bx])[:,np.newaxis]
        By = np.array([By])[:,np.newaxis]
        Bz = np.array([Bz])[:,np.newaxis]
        self.cosmu = vx.flatten()*Bx+vy.flatten()*By+vz.flatten()*Bz
        self.cosmu = self.cosmu.reshape(self.cosmu.shape[0],vdim[0],vdim[1])[0,:,:]   # removing unnecessary dimension of array through slicing
        self.mu = np.arccos(self.cosmu)/np.pi*180.

    def _calc_pitch_ang_cov(self,cosmu_b = np.arange(-1.,1.01,0.01),mu_b = np.arange(0,180.1,1.)):
        tdim = self.B.shape[1]
        self.cosmu_cov = np.zeros((tdim,16,cosmu_b.shape[0]-1))
        self.cosmu_cov_rel = np.zeros((tdim,16,cosmu_b.shape[0]-1))
        self.cosmu_cov_tot_rel = np.zeros((tdim,16,cosmu_b.shape[0]-1))
        self.mu_cov = np.zeros((tdim,16,mu_b.shape[0]-1))
        self.mu_cov_rel = np.zeros((tdim,16,mu_b.shape[0]-1))
        self.mu_cov_tot_rel = np.zeros((tdim,16,mu_b.shape[0]-1))
        for t in range(tdim):
            for i in range(16):
                x,y = np.histogram(self.cosmu[t].flatten(),cosmu_b,weights=self.geometry_factor[i].flatten()*self.scale.flatten())
                self.cosmu_cov[t,i] = x
                self.cosmu_cov_rel[t,i] = x/(self.geometry_factor[i].flatten()*self.scale.flatten()).sum()
                x,y = np.histogram(self.mu[t].flatten(),mu_b,weights=self.geometry_factor[i].flatten()*self.scale.flatten())
                self.mu_cov[t,i] = x
                self.mu_cov_rel[t,i] = x/(self.geometry_factor[i].flatten()*self.scale.flatten()).sum()
            for i in range(16):
                self.cosmu_cov_tot_rel[t,i] = self.cosmu_cov[t,i]/self.cosmu_cov[t].sum(axis=0) 
                self.mu_cov_tot_rel[t,i] = self.mu_cov[t,i]/self.mu_cov[t].sum(axis=0) 
        self.cosmu_cov_tot = self.cosmu_cov.sum(axis = 1)
        self.mu_cov_tot = self.mu_cov.sum(axis = 1)
    def _set_pixel_pos(self):
        """
        Here the central position of all pixels for all 16 pixel are defined in x,y in plane of pixels
        """
        # all Pixel x-pos (3D)
        self.pixel_pos[:,0] = self.d
        
        # Pixel 0 = Background pixel
        self.pixel_pos[0,1] = (4.+2*0.46+1.17)
        self.pixel_pos[0,2] = 0.

        # z-positions columnwise (2D - y-position)
        for i in range(1,6):
            self.pixel_pos[i::5,2] = (3.-i)*(2.+0.46+0.77)
        # y-position rowwise (2D - x-position)
        for i in range(0,3):
            self.pixel_pos[i*5+1:i*5+6,1] = (1-i)*(2.+0.46)
        
    def _set_pinhole(self):
        ph = np.meshgrid(np.linspace(-0.7,0.7,self.npy),np.linspace(-1.4,1.4,self.npz))
        self.pinhole = np.zeros((3,self.npz,self.npy))
        self.pinhole[0] = self.d
        self.pinhole[1] = ph[0]
        if self.mag:
            self.pinhole[1] += 0.25
        self.pinhole[2] = ph[1]

    def _calc_pA(self):
        """
        Calculates projected Area of the pinhole under all viewing directions
        """
        for i,phi in enumerate(self.phi_r):
            for j,theta in enumerate(self.theta_r):
                self.pA_pinhole[i,j] = np.cos(phi)*np.cos(theta)*self.A_pinhole
        
        
    def _calc_directions(self):
        for i,phi in enumerate(self.phi_r):
            for j,theta in enumerate(self.theta_r):
                y,z = self._calc_direction(phi,theta)
                self.directions[:,1,i,j] = y
                self.directions[:,2,i,j] = z
        self.directions[:,0,:,:] = self.d
        
    def _calc_direction(self,phi,theta):
        y = np.tan(phi)*self.d
        z = np.tan(theta)*np.sqrt((np.tan(phi)**2+1)*self.d**2)
        py = self.pinhole[1] + y
        pz = self.pinhole[2] + z
        return py.flatten(),pz.flatten()

    def _calc_geometry_factor(self):
        y = self.directions[:,1,:,:]
        z = self.directions[:,2,:,:]
        for p in range(16):
            py = self.pixel_pos[p,1]
            pz = self.pixel_pos[p,2]
            mask = (y>=py-1.)*(y<=py+1.)*(z>=pz-1.)*(z<=pz+1.)
            hitfrac = mask.sum(axis=0)/self.dim
            self.hitfrac = hitfrac
            self.geometry_factor[p] = hitfrac * self.pA_pinhole

class VDSTEP(object):
    def __init__(self, angstep = 1., B = np.array([np.cos(4/18*np.pi),-np.sin(4/18.*np.pi),0.0]), mag = False):
        """
        Coordinatesystem :
        x -> Axis through center of STEP (i.e. Pixel 8) orientation from pinhole to detector array (0,0,0 = center of Pinhole)
        y -> from center to background pixel (i.e. parallel to T in RTN)
        Z -> completes the right handed triade (i.e. parallel to N in RTN)
        """
        self.phi_d = np.arange(-16,16.1,angstep)    # angles in Instrument Coordinates
        self.theta_d = np.arange(-30.,30.1,angstep) # angles in Instrument Coordinates
                                                    # theta is not in spherical coordinates but against the main plane, i.e. centre of bckgrd and pixel 3,8,13
        self.phi_r = self.phi_d/180.*np.pi          # calculating angles in radians
        self.theta_r = self.theta_d/180.*np.pi      # calculating angles in radians
        self.phi_offs_d = -35.
        self.theta_offs_d = 0.
        self.phi_offs_r = self.phi_offs_d/180.*np.pi
        self.theta_offs_r = self.theta_offs_d/180.*np.pi
        self.d = 15 # Distance between pinhole(front) and plane of the ssd-array in mm
        self.B = B
        self.directions = np.zeros((3,self.phi_r.shape[0],self.theta_r.shape[0]))
        self.vdir = np.zeros((3,self.phi_r.shape[0],self.theta_r.shape[0]))
        self.cosmu = np.zeros((self.phi_r.shape[0],self.theta_r.shape[0]))
        self.mu = np.zeros((self.phi_r.shape[0],self.theta_r.shape[0]))
        self.A_pinhole = 0.14 * 0.28  # cm^2
        self.pA_pinhole = np.zeros((self.phi_r.shape[0],self.theta_r.shape[0]))
        self.hitfrac = np.zeros((16,self.phi_r.shape[0],self.theta_r.shape[0]))
        self.hitpoints = np.zeros((16,self.phi_r.shape[0],self.theta_r.shape[0],4,3))
        self.geometry_factor = np.zeros((16,self.phi_r.shape[0],self.theta_r.shape[0]))
        self.scale = np.ones((self.phi_r.shape[0],self.theta_r.shape[0]))*np.cos(self.theta_r)[np.newaxis,:]
        self.pixel_pos = np.zeros((16,3)) 
        self.mag = mag
        self._set_pixel_pos()
        self._set_pinhole()
        self._calc_pA()
        self._calc_directions()
        self._calc_geometry_factor()
        self._calc_vdir()
        self._calc_pitch_ang()
        #self._calc_pitch_ang_cov()
        
    def _calc_vdir(self):
        # calc unity velocity vectors in SUN-SC Coordinates
        for i,p in enumerate(self.phi_r):
            for j,t in enumerate(self.theta_r):
                th = -t + np.pi/2. + self.theta_offs_r
                ph = p + self.phi_offs_r
                self.vdir[0,i,j] = np.sin(th)*np.cos(ph)
                self.vdir[1,i,j] = np.sin(th)*np.sin(ph)
                self.vdir[2,i,j] = np.cos(th)

    def _calc_pitch_ang(self):
        Bx = self.B[0]/np.sqrt(self.B[0]**2+self.B[1]**2+self.B[2]**2)
        By = self.B[1]/np.sqrt(self.B[0]**2+self.B[1]**2+self.B[2]**2)
        Bz = self.B[2]/np.sqrt(self.B[0]**2+self.B[1]**2+self.B[2]**2)
        vx = self.vdir[0]
        vy = self.vdir[1]
        vz = self.vdir[2]
        vdim = vx.shape
        Bx = np.array([Bx])[:,np.newaxis]
        By = np.array([By])[:,np.newaxis]
        Bz = np.array([Bz])[:,np.newaxis]
        self.cosmu = vx.flatten()*Bx+vy.flatten()*By+vz.flatten()*Bz
        self.cosmu = self.cosmu.reshape(self.cosmu.shape[0],vdim[0],vdim[1])[0,:,:]   # removing unnecessary dimension of array through slicing
        self.mu = np.arccos(self.cosmu)/np.pi*180.

    def _calc_pitch_ang_cov(self,cosmu_b = np.arange(-1.,1.01,0.01),mu_b = np.arange(0,180.1,1.)):
        tdim = self.B.shape[1]
        self.cosmu_cov = np.zeros((tdim,16,cosmu_b.shape[0]-1))
        self.cosmu_cov_rel = np.zeros((tdim,16,cosmu_b.shape[0]-1))
        self.cosmu_cov_tot_rel = np.zeros((tdim,16,cosmu_b.shape[0]-1))
        self.mu_cov = np.zeros((tdim,16,mu_b.shape[0]-1))
        self.mu_cov_rel = np.zeros((tdim,16,mu_b.shape[0]-1))
        self.mu_cov_tot_rel = np.zeros((tdim,16,mu_b.shape[0]-1))
        for t in range(tdim):
            for i in range(16):
                x,y = np.histogram(self.cosmu[t].flatten(),cosmu_b,weights=self.geometry_factor[i].flatten()*self.scale.flatten())
                self.cosmu_cov[t,i] = x
                self.cosmu_cov_rel[t,i] = x/(self.geometry_factor[i].flatten()*self.scale.flatten()).sum()
                x,y = np.histogram(self.mu[t].flatten(),mu_b,weights=self.geometry_factor[i].flatten()*self.scale.flatten())
                self.mu_cov[t,i] = x
                self.mu_cov_rel[t,i] = x/(self.geometry_factor[i].flatten()*self.scale.flatten()).sum()
            for i in range(16):
                self.cosmu_cov_tot_rel[t,i] = self.cosmu_cov[t,i]/self.cosmu_cov[t].sum(axis=0) 
                self.mu_cov_tot_rel[t,i] = self.mu_cov[t,i]/self.mu_cov[t].sum(axis=0) 
        self.cosmu_cov_tot = self.cosmu_cov.sum(axis = 1)
        self.mu_cov_tot = self.mu_cov.sum(axis = 1)
    def _set_pixel_pos(self):
        """
        Here the central position of all pixels for all 16 pixel are defined in x,y in plane of pixels
        """
        # all Pixel x-pos (3D)
        self.pixel_pos[:,0] = self.d
        
        # Pixel 0 = Background pixel
        self.pixel_pos[0,1] = (4.+2*0.46+1.17)
        self.pixel_pos[0,2] = 0.

        # z-positions columnwise (2D - y-position)
        for i in range(1,6):
            self.pixel_pos[i::5,2] = (3.-i)*(2.+0.46+0.77)
        # y-position rowwise (2D - x-position)
        for i in range(0,3):
            self.pixel_pos[i*5+1:i*5+6,1] = (1-i)*(2.+0.46)
        
    def _set_pinhole(self):
        self.pinhole = np.zeros((3))
        if self.mag:
            self.pinhole[1] = 0.25
    def _calc_pA(self):
        """
        Calculates projected Area of the pinhole under all viewing directions
        """
        for i,phi in enumerate(self.phi_r):
            for j,theta in enumerate(self.theta_r):
                self.pA_pinhole[i,j] = np.cos(phi)*np.cos(theta)*self.A_pinhole
        
        
    def _calc_directions(self):
        for i,phi in enumerate(self.phi_r):
            for j,theta in enumerate(self.theta_r):
                y,z = self._calc_direction(phi,theta)
                self.directions[1,i,j] = y
                self.directions[2,i,j] = z
        self.directions[0,:,:] = self.d
        
    def _calc_direction(self,phi,theta):
        y = np.tan(phi)*self.d
        z = np.tan(theta)*np.sqrt((np.tan(phi)**2+1)*self.d**2)
        py = self.pinhole[1] + y
        pz = self.pinhole[2] + z
        return py,pz

    def _calc_geometry_factor(self):
        for i in range(self.geometry_factor.shape[1]):
            for j in range(self.geometry_factor.shape[2]):
                y = self.directions[1,i,j]
                z = self.directions[2,i,j]
                for p in range(16):
                    
                    dy = min(self.pixel_pos[p,1]+1,y+0.7) - max(self.pixel_pos[p,1]-1,y-0.7)
                    dz = min(self.pixel_pos[p,2]+1,z+1.4) - max(self.pixel_pos[p,2]-1,z-1.4)
                    if (dy>0) * (dz>0):
                        hitfrac = dy*dz /100. / self.A_pinhole
                        self.hitpoints[p,i,j,:,0] = self.d
                        self.hitpoints[p,i,j,0,1] = min(self.pixel_pos[p,1]+1,y+0.7)
                        self.hitpoints[p,i,j,0,2] = min(self.pixel_pos[p,2]+1,z+1.4)
                        self.hitpoints[p,i,j,1,1] = max(self.pixel_pos[p,1]-1,y-0.7)
                        self.hitpoints[p,i,j,1,2] = min(self.pixel_pos[p,2]+1,z+1.4)
                        self.hitpoints[p,i,j,2,1] = max(self.pixel_pos[p,1]-1,y-0.7)
                        self.hitpoints[p,i,j,2,2] = max(self.pixel_pos[p,2]-1,z-1.4)
                        self.hitpoints[p,i,j,3,1] = min(self.pixel_pos[p,1]+1,y+0.7)
                        self.hitpoints[p,i,j,3,2] = max(self.pixel_pos[p,2]-1,z-1.4)
                    else:
                        hitfrac = 0.
                    self.hitfrac[p,i,j] = hitfrac
                    self.geometry_factor[p,i,j] = hitfrac * self.pA_pinhole[i,j]

def plot_cosmu_cov(vd,t=0):
    fig,ax = create_FoV_frame()
    cmap = cm.jet
    CS = ax[0].contourf(vd.phi_d+35,vd.theta_d,vd.mu[t].T[-1::-1,-1::-1],np.arange(0.,180.01,1.),cmap=cmap)
    CS2 = ax[0].contour(vd.phi_d+35,vd.theta_d,vd.mu[t].T[-1::-1,-1::-1],np.arange(0.,180.01,5.),colors='black')
    ax[0].clabel(CS2, CS2.levels, inline=True, fontsize=10,fmt = '%3.0f')
    bphi= np.arctan2(vd.B[1][t],vd.B[0][t])/np.pi*180.
    btheta= np.arctan2(vd.B[2][t],np.sqrt(vd.B[0][t]**2+vd.B[1][t]**2))/np.pi*180.
    print(bphi,btheta)
    ax[0].plot(-bphi,-btheta,"x",color='k')
    cb_ax = fig.add_axes([0.1, 0.94, 0.8, 0.02])
    cb_ax.xaxis.set_label_position('top')
    cbar = fig.colorbar(CS, cax = cb_ax, orientation = 'horizontal')
    cbar.set_ticks(np.arange(0.,180.1,10.))
    cbar.set_label('Pitch Angle Whole FoV | Pixel FoV', labelpad = -45)
    for i in range(1,16):
        tmu = 1.*vd.mu[t]#/np.pi*180.
        tmu[vd.geometry_factor[i]==0]=np.NaN
        ax[i].contourf(vd.phi_d+35,vd.theta_d,tmu.T[-1::-1,-1::-1],np.arange(-0.,180.01,5.),cmap=cmap)
        CS2 = ax[i].contour(vd.phi_d+35,vd.theta_d,vd.mu[t].T[-1::-1,-1::-1],np.arange(0.,180.01,5.),colors="black",linestyles = 'dotted')
        ax[i].plot(-bphi,-btheta,"x",color='k')

def plot_pic_frac(vd):
    fig,ax = create_FoV_frame()

    ns = vd.geometry_factor.sum(axis=0)
    rc = np.zeros(vd.geometry_factor.shape)
    nc = vd.geometry_factor/ns
    cmap = cm.jet
    nc[nc==0.]=np.NaN
    ns = ns/np.amax(ns)
    ns[ns==0]=np.NaN
    ax[0].contourf(vd.phi_d+35,vd.theta_d,ns.T[-1::-1,-1::-1],np.arange(0,1.01,0.1),cmap=cmap)
    for i in range(1,16):
        CS = ax[i].contourf(vd.phi_d+35,vd.theta_d,nc[i].T[-1::-1,-1::-1],np.arange(0.,1.01,.1),cmap = cmap,vmin=0.,vmax=1.)
    cb_ax = fig.add_axes([0.1, 0.94, 0.8, 0.02])
    cb_ax.xaxis.set_label_position('top')
    cbar = fig.colorbar(CS, cax = cb_ax, orientation = 'horizontal')
    cbar.set_label('Normalised Geometryfactor | Pixel Contribution', labelpad = -45)
    return fig,ax
    
def plot_illumination(vd,save = False):
    fig,ax = plot_pic_frac(vd)
    for i in range(16):
        ax[i].plot([50,33],[12,29],'x-',color='k')
    if save:
        plt.savefig('/home/ivar/berger/projects/solo/plots/Illumination.png')
        plt.savefig('/home/ivar/berger/projects/solo/plots/Illumination.pdf')
    
def plot_FoV(vd,phii = 0,thetai = 0, plist = range(16)):
    '''phii and thetai are the indices of the arrays of the angles.'''
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_zlim(-10.,10.)
    ax.set_ylim(-10.,10.)
    ax.set_xlim(-5.,15)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.view_init(0,180)

    phy = np.array([-0.7,0.7,0.7,-.7,-.7]) + vd.pinhole[1]
    phz = np.array([-1.4,-1.4,1.4,1.4,-1.4]) + vd.pinhole[0]
    phx = np.array([0.,0.,0.,0.,0.]) + vd.pinhole[2]
    ax.plot(phx,phy,phz,color="k")
    for p in plist:
        x = vd.pixel_pos[p,0]
        y = vd.pixel_pos[p,1]
        z = vd.pixel_pos[p,2]
        py = np.array([y-1.,y+1,y+1,y-1,y-1])
        pz = np.array([z-1.,z-1.,z+1,z+1,z-1])
        px = np.array([x,x,x,x,x])
        ax.plot(px,py,pz,color="k")
        ax.text(x,y,z,'%i'%(p),ha = 'center', va = 'center')
    y = vd.directions[1,phii,thetai]
    z = vd.directions[2,phii,thetai]
    x = vd.directions[0,phii,thetai]
    y0 = phy
    z0 = phz
    x0 = phx
    y = y0 + y - vd.pinhole[1] 
    z = z0 + z - vd.pinhole[0] 
    x = x0 + x - vd.pinhole[2] 
    dy = (y-y0)/15.
    dz = (z-z0)/15.
    dx = -5*np.ones((x.shape))
    y5 = y0 + dy * dx
    z5 = z0 + dz * dx
    for i in range(x.shape[0]-1):
        ax.plot([x0[i],x[i]],[y0[i],y[i]],[z0[i],z[i]],color="r",alpha=0.5)
        ax.plot([x0[i],dx[i]],[y0[i],y5[i]],[z0[i],z5[i]],color="r",alpha=0.5)
        ax.plot([x[i],x[i+1]],[y[i],y[i+1]],[z[i],z[i+1]],color="r",alpha=1.)
    #verts = [list(zip(x[:-1],y[:-1],z[:-1]))]
    #ax.add_collection3d(Poly3DCollection(verts,facecolor='r',alpha = 0.5))   
    for p in plist:
        if sum(vd.hitpoints[p,phii,thetai,:,0]):
            #x = np.append(vd.hitpoints[p,phii,thetai,:,0],vd.hitpoints[p,phii,thetai,0,0])
            #y = np.append(vd.hitpoints[p,phii,thetai,:,1],vd.hitpoints[p,phii,thetai,0,1])
            #z = np.append(vd.hitpoints[p,phii,thetai,:,2],vd.hitpoints[p,phii,thetai,0,2])
            x = vd.hitpoints[p,phii,thetai,:,0]
            y = vd.hitpoints[p,phii,thetai,:,1]
            z = vd.hitpoints[p,phii,thetai,:,2]
            verts = [list(zip(x,y,z))]
            #ax.plot(x,y,z,color="g",alpha=1.)
            ax.add_collection3d(Poly3DCollection(verts,facecolor='g'))   
    
    #ax.plot([np.amin(x),np.amax(x),np.amax(x),np.amin(x),np.amin(x)],[np.amin(y),np.amin(y),np.amax(y),np.amax(y),np.amin(y)],[15,15,15,15,15],color="r")

    #x = vd.directions[:,0,int(vd.phi_r.shape[0]/2),int(vd.theta_r.shape[0]/2)]
    #y = vd.directions[:,1,int(vd.phi_r.shape[0]/2),int(vd.theta_r.shape[0]/2)]
    #z = 15. * ones((x.shape))
    #x0 = vd.pinholex.flatten()
    #y0 = vd.pinholey.flatten()
    #z0 = zeros((x.shape))
    #dx = (x-x0)/15.
    #dy = (y-y0)/15.
    #dz = -5*ones((x.shape))
    #x5 = x0 + dx * dz
    #y5 = y0 + dy * dz
    #for i in range(x.shape[0]):
    #    ax.plot([x0[i],x[i]],[y0[i],y[i]],[z0[i],z[i]],color="g",alpha=0.1)
    #    ax.plot([x0[i],x5[i]],[y0[i],y5[i]],[z0[i],dz[i]],color="g",alpha=0.1)
    #ax.plot([amin(x),amax(x),amax(x),amin(x),amin(x)],[amin(y),amin(y),amax(y),amax(y),amin(y)],[15,15,15,15,15],color="g")
    return ax


def plot_geomf(vd):
    phib = vd.phi_d-35
    #phib = phib[-1::-1]
    thetab = vd.theta_d
    fig = plt.figure()
    for j in range(3):
        for i in range(5):
            fig.add_subplot(5,3,1+j+i*3)
            plt.contour(phib,thetab,vd.geometry_factor[i+1+j*5].T)
