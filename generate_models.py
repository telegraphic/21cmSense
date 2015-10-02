import numpy as np
import pylab as plt

def closest(arr, val):
     return np.argmin(np.abs(arr - val))

z   = np.genfromtxt('anastasia/z.txt')
k   = np.genfromtxt('anastasia/K.txt')
psh = np.genfromtxt('anastasia/PSh.txt')
pss = np.genfromtxt('anastasia/PSs.txt')

plt.imshow(psh, extent=(k[0], k[-1], z[0], z[-1]), aspect='auto', cmap='jet')
plt.colorbar(label='mK$^2$')
plt.xlabel('k (1/Mpc)')
plt.ylabel('Redshift (z)')
plt.savefig('anastasia-model.pdf')
plt.show()

z_targets = np.arange(10, 40, 0.25)

for z_target in z_targets:
    z_idx    = closest(z, z_target)

    model_out = np.column_stack((k, psh[z_idx], np.zeros_like(k)))
    np.savetxt('models/fialkov_z%2.2f_psh.txt' % z_target, model_out, delimiter='\t')

    plt.figure('Model k vs mK')
    plt.plot(model_out[:, 0], model_out[:, 1], c='#333333')
    plt.xlabel('k [1/Mpc]')
    plt.ylabel('mK$^2$')
    plt.title('z = 20.5')
    plt.savefig('figures/fialkov_z%2.2f_psh.pdf' % z_target)
    plt.clf()


