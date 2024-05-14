
import numpy as np
import h5py

def read_matlab_data():

    # Read .mat file
    f = h5py.File('BEM_data.mat','r')
    # print(f.keys())

    # Load the BEM variables that you need
    BEM_rR = np.asarray(f.get('rR'))
    BEM_alpha = np.asarray(f.get('alpha'))
    BEM_phi = np.asarray(f.get('phi'))
    BEM_rho = np.asarray(f.get('rho'))
    BEM_Ax = np.asarray(f.get('Ax'))
    BEM_Az = np.asarray(f.get('Az'))
    BEM_Gamma = np.asarray(f.get('gamma'))
    BEM_CT = np.asarray(f.get('CT'))
    BEM_CP = np.asarray(f.get('CP'))
    BEM_a = np.asarray(f.get('a'))
    BEM_aline = np.asarray(f.get('aprime'))
    BEM_vinf = np.asarray(f.get('uInf'))
    BEM_radius = np.asarray(f.get('rRotor'))

    # Close .mat file
    f.close()

    return [BEM_rR, BEM_alpha, BEM_phi, BEM_rho, BEM_Ax, BEM_Az, BEM_Gamma , BEM_CT, BEM_CP, BEM_a, BEM_aline, BEM_vinf, BEM_radius]

# Polar plot

# minVal = min([solverS0.alpha(:); solverS15.alpha(:); solverS30.alpha(:)]);
# maxVal = max([solverS0.alpha(:); solverS15.alpha(:); solverS30.alpha(:)]);
# t1 = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
# nexttile
# pplot = contourf(x, y, rad2deg(solverS0.alpha));
# h=colorbar;
# caxis manual
# caxis([rad2deg(minVal), rad2deg(maxVal)])
# xlabel(h,'\alpha','Rotation',0,'FontSize',18)
# xlabel('x/R (-)', 'FontSize',17)
# ylabel('y/R (-)', 'FontSize',17)
# colormap('default')
