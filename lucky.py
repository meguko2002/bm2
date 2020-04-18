import matplotlib.pyplot as plt
import numpy as np


class Sim:
    def __init__(self):
        self.t = np.arange ( 0.0 , 2.0 , 0.01 )
        self.z = np.cos ( 2 * np.pi * self.t )
        self.v = np.sin ( 2 * np.pi * self.t )
        self.fig , (self.ax1 , self.ax2) = plt.subplots ( nrows=2 , sharex=True )
        self.set_ax()
        self.ax1.set_xlim ( 0 , 2 )
        self.ax1.set_ylim ( -2 , 2 )
        self.ax2.set_xlim ( 0 , 2 )
        self.ax2.set_ylim ( -2 , 2 )


    def set_ax(self):
        self.ln_1 , = self.ax1.plot ( [] , [] , 'green' )
        self.ln_1.set_data ( self.t , self.v )
        self.ln_2 , = self.ax1.plot ( [] , [] , 'black' )
        self.ln_2.set_data ( self.t , self.v )

        self.ln_3 , = self.ax2.plot ( [] , [] , 'green' )
        self.ln_3.set_data (self.t, self.z)
        self.ln_4 , = self.ax2.plot ( [] , [] , 'black' )
        self.ln_4.set_data ( self.t , self.v )

    def onclick(self, event):
        x, y = event.xdata, event.ydata
        if event.button == 3:
            u = (x - 1) * self.v
            self.ln_2.set_data ( self.t , u)
            self.ln_4.set_data( self.t, self.v - y)
        elif event.button == 1:
            self.ln_2.remove()
            self.ln_4.remove()
            self.set_ax()

        self.fig.canvas.draw()

    def run(self):
        self.fig.canvas.mpl_connect ( 'button_press_event' , self.onclick )
        plt.show ()


sim = Sim ()
sim.run ()
