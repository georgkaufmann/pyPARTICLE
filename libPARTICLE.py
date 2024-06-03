"""
pyPARTICLE
library for particle transport model
2024-04-21
Georg Kaufmann
"""

import numpy as np
import matplotlib.pyplot as plt
import sys,os
import libPARTICLE


#================================#
def particleInit(Nfloat=50,Nsettled=10,**kwargs):
    """
    Particle motion
    input:
      Nfloat,Nsettled    - number of floating/settled particles
    **kwargs:
      sidex,sidey        - x- and y-dimension of box
      seedx,seedy        - center of floating particle cloud
      sigmax,sigmay      - standard deviation for random location
    output:
      pX,pY              - x- and y-coordinates of particles
      pState             - state of particles (0-settled,>0-floating)
      pMaterial          - material of particle
    """
    # default settings
    sidex  = 1.0; sidey  = 0.5
    seedx  = 0.5; seedy  = 0.4
    sigmax = 0.05; sigmay = 0.05
    # test for user settings
    print('User settings: ',kwargs)
    for i in kwargs:
        if (i=='sidex'): sidex = kwargs[i]
        if (i=='sidey'): sidey = kwargs[i]
        if (i=='seedx'): seedx = kwargs[i]
        if (i=='seedy'): seedy = kwargs[i]
        if (i=='sigmax'): sigmax = kwargs[i]
        if (i=='sigmay'): sigmay = kwargs[i]
    # total number of particles
    N = Nfloat + Nsettled
    # random location for floating particles
    rng    = np.random.default_rng(seed=12)
    pX     = rng.uniform(seedx-sigmax,seedx+sigmax,Nfloat)
    pY     = rng.uniform(seedy-sigmay,seedy+sigmay,Nfloat)
    pState    = np.ones(Nfloat,dtype='int')
    pMaterial = np.ones(Nfloat,dtype='int')
    # set coordinates of settled particles
    pX = np.append(pX,np.linspace(0,sidex,Nsettled),axis=0)
    pY = np.append(pY,np.zeros(Nsettled),axis=0)
    pState    = np.append(pState,np.zeros(Nsettled,dtype='int'),axis=0)
    pMaterial = np.append(pMaterial,np.zeros(Nsettled,dtype='int'),axis=0)
    return pX,pY,pState,pMaterial
    

#================================#
def particlePlot(pX,pY,pState,pMaterial,time=0,isaved=0,show=False,**kwargs):
    """
    Particle motion
    input:
      pX,pY              - x- and y-coordinates of particles
      pState             - state of particles (0-settled,>0-floating)
      pMaterial          - material of particle
    **kwargs:
      sidex,sidey        - x- and y-dimension of box
      path               - path for figure and data files
      name               - name of model
    output:
      (to figure)
    """
    colors = ['black','red','green','blue']
    # default settings
    sidex  = 1.0; sidey  = 0.5
    path  = 'div/test1'
    name  = 'test1'
    # test for user settings
    for i in kwargs:
        if (i=='sidex'): sidex = kwargs[i]
        if (i=='sidey'): sidey = kwargs[i]
        if (i=='path'):  path = kwargs[i]
        if (i=='name'):  name = kwargs[i]
    # check for directory for plotting
    if not os.path.isdir(path):
        os.mkdir(path)
    # plot
    offsetx = 0.1*sidex; offsety = 0.1*sidey
    filename = name+f"-{isaved:04}.png"
    
    fig,ax = plt.subplots(1,1,figsize=(6,6))
    ax.set_title('time: '+str(round(time,2)))
    ax.set_xlim([0-offsetx,sidex+offsetx])
    ax.set_ylim([0-offsety,sidey+offsety])
    ax.set_aspect('equal', 'box')
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    for i in range(pX.shape[0]):
        if (pState[i]==0):
            plt.plot(pX[i],pY[i],lw=0,marker='*',markersize=4,color=colors[pMaterial[i]],alpha=1.0)
        if (pState[i]>0):
            plt.plot(pX[i],pY[i],lw=0,marker='o',markersize=4,color=colors[pMaterial[i]],alpha=0.4)
    plt.savefig(path+'/'+filename)
    if (not show):
        plt.close()
    return


#================================#
def particleMaterials(show=False):
    """
    Particle motion
    Define material properties
    input:
      (none)
    output:
      materials - json dictionary for materials
    """
    materials = [
    {"id": 0, "type": 'settled', "radius":    1e-3, "density": 2600},
    {"id": 1, "type": 'gravel',  "radius":    1e-3, "density": 2600},
    {"id": 2, "type": 'sand',    "radius":  0.4e-3, "density": 2400},
    {"id": 3, "type": 'silt',    "radius": 0.10e-3, "density": 2200}
    ]
    if (show):
        print('Number of materials: {}'.format(len(materials)))
        for material in materials:
            print("id: {}\ntype: {}\nradius: {} m\ndensity: {} kg/3\n".format(material['id'],material['type'],material['radius'],material['density']))
    return materials


#================================#
def particleRun(pX,pY,pState,pMaterial,materials,tmax=10,dt=1,TC=10,**kwargs):
    """
    Particle motion
    """
    # default settings
    sidex = 1.0; sidey  = 0.5
    xdiff = 0.01; ydiff = 0.01
    path  = 'div/test1'
    name  = 'test1'
    # test for user settings
    print('User settings: ',kwargs)
    for i in kwargs:
        if (i=='sidex'): sidex = kwargs[i]
        if (i=='sidey'): sidey = kwargs[i]
        if (i=='path'): path = kwargs[i]
        if (i=='name'): name = kwargs[i]
        if (i=='xdiff'): xdiff = kwargs[i]
        if (i=='ydiff'): ydiff = kwargs[i]
        if (i=='xadv'):  xadv = kwargs[i]
        if (i=='yadv'):  yadv = kwargs[i]
        if (i=='tmax'): tmax = kwargs[i]
        if (i=='dt'): dt = kwargs[i]
        if (i=='TC'): TC = kwargs[i]
    # check for directory for plotting
    if not os.path.isdir(path):
        os.mkdir(path)
    print('path: ',path)
    print('name: ',name)
    # mean radius of all particles
    meanRadius = 0.
    for i in range(len(materials)):
        meanRadius += materials[i]['radius']
    meanRadius = meanRadius / len(materials)
    print('mean radius {} m'.format(meanRadius))
    # initialize random-number generator
    rng = np.random.default_rng(seed=12)
    t   = dt
    isaved = 1
    while (t < tmax):
        print(isaved,end=' ')
        # mark floating and settled particles in two different arrays
        pFloating = np.array([],dtype='int')
        pSettling = np.array([],dtype='int')
        # loop over particles
        for i in range(pX.shape[0]):
            if (pState[i] > 0):
                # advect particle
                radius = materials[pMaterial[i]]['radius']
                pX[i] = pX[i] + libPARTICLE.particleVelocityShear(sidey-pY[i],sidey,v0=xadv,TC=TC)*dt
                if (yadv):    
                    pY[i] = pY[i] - libPARTICLE.particleVelocityStokes(radius,TC=TC)*dt
                # diffuse particle
                #pX[i] = pX[i] + xdiff*rng.uniform(-1,1)
                #pY[i] = pY[i] + ydiff*rng.uniform(-1,1)
                pX[i] = pX[i] + meanRadius*np.random.normal(loc=0,scale=1)
                pY[i] = pY[i] + meanRadius*np.random.normal(loc=0,scale=1)
            if (pState[i]==1): pFloating = np.append(pFloating,[i])
            if (pState[i]==0): pSettling = np.append(pSettling,[i])
        # keep fluid particles within domain
        pX = np.where(pX < 0, 0, pX)
        pX = np.where(pX > sidex, sidex, pX)
        pY = np.where(pY < 0, 0, pY)
        pY = np.where(pY > sidey, sidey, pY)
        # check, if floating particle is close to settled particle, if yes, attach
        for i,ip in enumerate(pFloating):
            for j,ic in enumerate(pSettling):
                dist = materials[pMaterial[ip]]['radius'] + materials[pMaterial[ic]]['radius']
                #dist = 0.01
                if (np.sqrt((pX[ip]-pX[ic])**2 + (pY[ip]-pY[ic])**2) < 2*dist):
                    ##print(i,ip,j,ic)
                    pState[ip] = 0
        # plot particle locations to file
        libPARTICLE.particlePlot(pX,pY,pState,pMaterial,time=t,isaved=isaved,sidex=sidex,sidey=sidey,path=path,name=name)
        t      += dt
        isaved += 1
    return


#================================#
def particleVelocityShear(z,h,dpdx=0.,v0=0,TC=10):
    """
    Particle motion: 
    Couette and/or Poiseuille shear-flow profile
    input:
      z        - vertical coordninate
      h        - fluid layer thickness
      dpdx     - pressure gradient
      v0       - horizontal velocity along surface 
    output:
      vShear   - horizontal shear velocity profile
    """
    vShear = 1/2/libPARTICLE.water_viscosity(TC) * dpdx * ( z**2 - z*h ) + v0*(1. - z/h)
    return vShear


#================================#
def particleVelocityStokes(r,rhoParticle=1650,TC=10):
    """
    Particle motion: 
    Stokes drag velocity
    """
    g = 9.81
    vStokes = 2*(rhoParticle - libPARTICLE.water_density(TC)) *g / 9. / libPARTICLE.water_viscosity(TC) * r**2
    return vStokes


#================================#
def water_viscosity(TC=20.):
    import numpy as np
    """
    ! water viscosity as a function of temperature
    ! Kestin, J., Sokolov, M., and Wakeham, W.: 
    ! Viscosity of Liquid Water in the Range -8C to 150C 
    ! J. Phys. Chem. Ref. Data, 7 (3), 941-948, 1978
    ! input:
    !  temperature     - [C]
    ! output:
    !  water_viscosity - [Pa s]
    """
    visc20      = 1002.e-6
    t20         = 20. - TC
    water_viscosity   = (20.-TC) / (96.+TC)* (1.2378 - 1.303e-3*t20 + 3.06e-6*t20**2 + 2.55e-8*t20**3)
    water_viscosity   = 10.**water_viscosity
    water_viscosity   = water_viscosity * visc20
    return water_viscosity


#================================#
def water_density(TC=20):
    import numpy as np
    """
    ! water density as a function of temperature
    ! Kell, G. S., 
    ! Effects of isotopic composition, temperature, 
    ! pressure, and dissolved gases on the density of liquid water 
    ! J. Phys. Chem. Ref. Data, 6, 1109, 1977.
    ! input:
    !  temperature     - [C]
    ! output:
    !  water_density   - [kg/m3]
    """
    water_density = (999.8437 + 67.8782e-3*TC + 103.1412e-6*TC**3 + 15.95835e-9*TC**5 + 636.8907e-15*TC**7) \
    / (1. + 9.090169e-6*TC**2 + 1.4511976e-9*TC**4 + 134.84863e-15*TC**6 + 2.008615e-18*TC**8)
    return water_density
#================================#