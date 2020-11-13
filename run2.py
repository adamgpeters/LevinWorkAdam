import numpy as np

# Neuron
from neuron import h

# visualization
import matplotlib.pyplot as plt

def plots(data1, data2):
    # Dividing figure
    x, y = [7, 5]
    fig, ax = plt.subplots(x, y, sharex=True, sharey=True, figsize=(50.0, 30.0))

    t = np.linspace(0, data1.shape[1] - 1, data1.shape[1])

    # Drawing images
    current = 20
    current_add = 10
    current_counter = 0
    for xi in range(x):
        for yi in range(y):
            # ax[xi, yi].title.set_text(f'Injection current {current + current_add * current_counter}pA')
            ax[xi, yi].plot(t, data1[current_counter, :], lw=2, label='observation')
            ax[xi, yi].plot(t, data2[current_counter, :], '--', lw=2, label='posterior sample')

            if yi == 0:
                ax[xi, yi].set(ylabel='mV')
            if xi == x - 1:
                ax[xi, yi].set(xlabel='Millisecond')

            current_counter += 1

    plt.subplots_adjust(hspace=0.3)

    plt.savefig('volts.png')
    plt.draw()
    plt.pause(0.5)

    return [fig, ax]
    
    
    
# Observed data
h.load_file("stdrun.hoc")  # for run control
h.load_file("neuron.hoc")  # run the model
h.initWT()
h.runWTsim()
h.initMT()
h.runMuTsim()
num_n = h.data_vecs_WT.__len__()
if num_n > 0:
    n_shape = [num_n, h.data_vecs_WT[0].__len__()]
    WT_voltage = np.zeros(n_shape, dtype=np.float)
    MT_voltage = np.zeros(n_shape, dtype=np.float)
    for y in range(n_shape[0]):
        WT_voltage[y, :] = h.data_vecs_WT[y].as_numpy()
        MT_voltage[y, :] = h.data_vecs_Mut[y].as_numpy()

plots(WT_voltage, MT_voltage)