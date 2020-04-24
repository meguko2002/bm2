import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons

fig, ax = plt.subplots()
plt.subplots_adjust(left=0.25, bottom=0.25)
t = np.arange(0.0, 5, 0.001)
f0 = 3
delta_f = 1
s = np.sin(2 * np.pi * f0 * t)
l, = plt.plot(t, s, lw=2)
ax.margins(x=0)

axcolor = 'lightgoldenrodyellow'
axfreq = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
sfreq = Slider(axfreq, 'Freq', 1, 10, valinit=f0, valstep=delta_f)


def update(val):
    # freq = sfreq.val
    # l.set_ydata(np.sin(2*np.pi*freq*t))
    ax.set_xlim(0, sfreq.val)
    fig.canvas.draw_idle()


sfreq.on_changed(update)

# resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
# button = Button(resetax, 'Reset' , color=axcolor, hovercolor='0.975')
#
# def reset(event):
#     sfreq.reset()
# button.on_clicked(reset)


plt.show()