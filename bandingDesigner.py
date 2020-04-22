import cv2
import matplotlib.patches as patches
import matplotlib.pylab as plt
from matplotlib.widgets import RadioButtons , RectangleSelector
import numpy as np
import imageanalyzer2
import PySimpleGUI as sg


class Simuration ():
    def __init__(self , file , dt , isimage=True):
        self.file = file
        self.dt = dt
        # self.ps = 226
        self.f_range = 100
        self.amp_max = 5
        self.isimage = isimage
        self.ref_timeobj = imageanalyzer2.Time ()
        self.sim_timeobj = imageanalyzer2.Time ()
        self.ref_fqobj = imageanalyzer2.Freq ()
        self.sim_fqobj = imageanalyzer2.Freq ()
        self.pre_fft = None
        self.ref_wave = None
        self.t_sq = self.ref_timeobj.t_sq
        self.fq_sq = self.ref_fqobj.fq_sq
        self.x_a , self.y_a , self.x_b , self.y_b = None , None , None , None
        self.filter = None
        self.selected_area = None  # or 'fft'
        self.set_fig ()
        self.set_ax ()
        self.set_frame ()

    def set_fig(self):
        if self.isimage:
            self.fig , (self.ax_original_data , self.ax_ref_image , self.ax_sim_image , \
                        self.ax_wave , self.ax_fft) = plt.subplots ( nrows=5 , figsize=(8 , 8) )
        else:
            self.fig , (self.ax_original_data , self.ax_wave , self.ax_ref_image , self.ax_sim_image , \
                        self.ax_fft) = plt.subplots ( nrows=5 , figsize=(8 , 8) )

    def set_ax(self):
        self.line_refwave , = self.ax_wave.plot ( [] , [] , 'green' , \
                                                  linestyle="dashed" , linewidth=0.5 , label='ref' )
        self.line_reffft , = self.ax_fft.plot ( [] , [] , 'green' , \
                                                linestyle="dashed" , linewidth=0.5 , label='ref' )
        self.line_simwave , = self.ax_wave.plot ( [] , [] , 'black' , \
                                                  linewidth=1 , label='sim' )
        self.line_simfft , = self.ax_fft.plot ( [] , [] , 'black' , \
                                                linewidth=1 , label='sim' )

    def set_frame(self):
        self.ax_original_data.set_title ( 'original image' )
        self.ax_ref_image.set_title ( 'ref image' )
        for spine in ('left' , 'right' , 'top' , 'bottom'):
            self.ax_ref_image.spines[spine].set_color ( "g" )
            self.ax_ref_image.spines[spine].set_linewidth ( 1 )
            self.ax_ref_image.spines[spine].set_linestyle ( "dashed" )
            self.ax_sim_image.spines[spine].set_color ( "black" )
            self.ax_sim_image.spines[spine].set_linewidth ( 2 )
        self.ax_sim_image.set_title ( 'sim image' )
        self.ax_wave.grid ( True )
        self.ax_fft.grid ( True )
        self.ax_wave.set_xlabel ( 'sec' )
        self.ax_wave.set_ylabel ( 'L*' )
        self.ax_fft.set_xlabel ( 'Hz' )
        self.ax_fft.set_ylabel ( 'amp' )
        self.ax_wave.legend ()
        self.ax_fft.legend ()

    def onselect(self , eclick , erelease):
        pass  # 座標を返すこともできるが、on,offの座標が逆転するバグがあるので使わない

    def onclick(self , event):
        if event.inaxes == self.ax_original_data:
            if event.button == 1:
                self.x_a , self.y_a = event.xdata , event.ydata
                self.selected_area = 'origin'
            elif event.button == 2 or event.button == 3:
                self.patch.set_visible ( not self.patch.get_visible () )
                self.fig.canvas.draw ()

        elif event.inaxes == self.ax_fft:
            if event.button == 1:
                self.x_a , self.y_a = event.xdata , event.ydata
                self.selected_area = 'fft'
            elif event.button == 2:  # reset
                self.sim_timeobj.set_wavedt ( self.ref_wave , dt )
                self.sim_fqobj.set_f ( self.ref_timeobj.fft )
                self.sim_draw ()
                self.fig.canvas.draw ()

            elif event.button == 3:  # undo
                tmp_fft = self.sim_fqobj.fft
                self.sim_fqobj.set_f ( self.pre_fft )
                self.pre_fft = tmp_fft
                self.sim_timeobj.set_wavedt ( self.sim_fqobj.inv_wave () , dt )
                self.sim_draw ()
                self.fig.canvas.draw ()

    def offclick(self , event):
        self.x_b , self.y_b = event.xdata , event.ydata
        if event.button == 1:
            if self.selected_area == 'origin' and event.inaxes == self.ax_original_data:
                self.ref_setting ()
            elif self.selected_area == 'fft' and event.inaxes == self.ax_fft:
                self.fft_change ()
            self.selected_area = None
            self.fig.canvas.draw ()

    def ref_setting(self):
        self.ax_original_data.cla ()
        if self.isimage:    self.set_original_image()
        else:   self.set_original_wave()
        self.trim_original_data()
        self.set_initial_wave ()
        self.ref_draw ()
        self.sim_draw ()
        self.set_initial_axis ()

    def trim_original_data(self):
        top = int ( min ( self.x_a , self.x_b ) )
        end = int ( max ( self.x_a , self.x_b ) )
        left = int ( min ( self.y_a , self.y_b ) )
        right = int ( max ( self.y_a , self.y_b ) )
        self.patch = patches.Rectangle ( xy=(top , left) , width=(end - top) , height=(right - left) , fill=False )
        self.ax_original_data.add_patch ( self.patch )
        if self.isimage:
            trimmed_image = self.raw_image[top: end , left: right]
            self.ref_wave = np.mean ( trimmed_image , axis=1 )
        else:
            self.ref_wave = self.original_wave[top: end]

    def set_initial_wave(self):
        self.ref_timeobj.set_wavedt ( self.ref_wave , dt )
        self.ref_fqobj.set_fdf ( self.ref_timeobj.fft , self.ref_timeobj.df )
        self.sim_fqobj = self.ref_fqobj
        self.sim_timeobj.set_wavedt ( self.sim_fqobj.inv_wave () , dt )
        self.t_sq = self.ref_timeobj.t_sq
        self.x_data = self.t_sq
        self.fq_sq = self.ref_fqobj.fq_sq

    def set_initial_axis(self):
        self.ax_wave.set_xlim ( 0 , self.ref_timeobj.n * self.ref_timeobj.dt )
        self.ax_wave.set_ylim ( 0 , 255 )
        self.ax_fft.set_xlim ( 0 , self.f_range )
        self.ax_fft.set_ylim ( 0 , self.amp_max )

    def ref_draw(self):
        self.line_refwave.set_data ( self.x_data , self.ref_timeobj.wave )
        self.line_reffft.set_data ( self.fq_sq , self.ref_fqobj.fft_abs )
        self.ax_ref_image.imshow ( self.ref_timeobj.image ( 1 ) , 'gray' , vmin=0 , vmax=255 , aspect='auto' )

    def sim_draw(self):
        self.line_simfft.set_data ( self.fq_sq , self.sim_fqobj.fft_abs )
        self.line_simwave.set_data ( self.x_data , self.sim_timeobj.wave )
        self.ax_sim_image.imshow ( self.sim_timeobj.image ( 1 ) , 'gray' , vmin=0 , vmax=255 , aspect='auto' )

    def fft_change(self):
        self.pre_fft = self.sim_fqobj.fft
        self.fft_filter ()
        self.sim_fqobj.set_f ( self.pre_fft * self.filter )
        self.sim_timeobj.set_wavedt ( self.sim_fqobj.inv_wave () , dt )
        self.line_simwave.set_data ( self.x_data , self.sim_timeobj.wave )
        self.line_simfft.set_data ( self.fq_sq , self.sim_fqobj.fft_abs )
        self.ax_sim_image.imshow ( self.sim_timeobj.image ( 1 ) , 'gray' , vmin=0 , vmax=255 , aspect='auto' )

    def fft_filter(self):
        prop = self.y_b / self.y_a
        out_freq = max ( self.x_a , self.x_b )
        in_freq = min ( self.x_a , self.x_b )
        self.filter = np.ones ( self.ref_timeobj.n )
        if in_freq > out_freq:
            tmp = in_freq
            in_freq = out_freq
            out_freq = tmp
        self.filter[(self.fq_sq > in_freq) & (self.fq_sq < out_freq)] = prop
        self.filter[(self.fq_sq > self.fq_sq[-1] - out_freq) & \
                    (self.fq_sq < self.fq_sq[-1] - in_freq)] = prop

    def set_original_image(self):
        self.raw_image = cv2.imdecode ( np.fromfile ( self.file , np.uint8 ) , 0 )
        self.ax_original_data.imshow ( self.raw_image.T , 'gray' , vmin=0 , vmax=255 , aspect='auto' )

    def set_original_wave(self):
        self.original_wave = np.loadtxt ( 'sinwave.txt' , dtype='uint8' )
        self.org_timeobj = imageanalyzer2.Time ( wave=self.original_wave , dt=self.dt )
        dot_sq = np.arange(0, self.org_timeobj.n)
        self.line_origilal_wave , = self.ax_original_data.plot (dot_sq, self.org_timeobj.wave, 'black' , \
                                                                 linewidth=1 , label='origin' )

    def run(self):
        if self.isimage:
            self.ax_original_data.set_title ( 'original image' )
            self.set_original_image ()
        else:
            self.ax_original_data.set_title ( 'original wave' )
            self.set_original_wave ()

        RS_im = RectangleSelector ( self.ax_original_data , self.onselect , drawtype='box' , useblit=True )
        RS_fft = RectangleSelector ( self.ax_fft , self.onselect , drawtype='box' , useblit=True )
        self.fig.canvas.mpl_connect ( 'button_press_event' , self.onclick )
        self.fig.canvas.mpl_connect ( 'button_release_event' , self.offclick )

        plt.show ()


if __name__ == '__main__':
    dt = 0.0005
    sim = Simuration ( ('ref.tif') , dt , isimage=0)
    sim.run ()
